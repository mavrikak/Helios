#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
HELIOS GUI — one window to run the full workflow

Covers:
- prepare / solve / post (wrapping run_sie.py)
- make-points (plane grids)
- cross-sections & near-field visualization (wrapping visualization.py with embedded Matplotlib)
- mesh viewer (mesh.mesh + optional interface planes from config.txt)

Requirements
------------
python -m pip install PySide6 matplotlib numpy

Usage
-----
python helios_gui.py --root <project_root>  # default: current folder

Notes
-----
- The GUI calls the existing project scripts via subprocess, preserving their
  behavior. Adjust the "Paths" in the Settings tab if your layout differs.
- Long-running tasks execute in a worker thread; stdout/stderr stream into the
  log pane. Outputs still land under sim_res/<sim>/ as usual.
"""
from __future__ import annotations
import sys, os, shlex, subprocess, threading, signal, time, re
from pathlib import Path
from dataclasses import dataclass

from PySide6 import QtCore, QtGui, QtWidgets
from PySide6.QtCore import Qt
from PySide6.QtGui import QDesktopServices
from PySide6.QtCore import QUrl

# Matplotlib for embedded plots
import matplotlib
matplotlib.use("QtAgg")
from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.image as mpimg
from matplotlib.collections import LineCollection
from matplotlib.patches import Patch

# reuse mesh helpers from pytools.plot_mesh
from pytools.plot_mesh import (
    read_mesh_mesh, auto_find_mesh, resolve_interfaces,
    build_transition_groups, decimate_indices, color_wheel, plot_mesh_3d
)

import numpy as np
import argparse

# ----------------------------- small utils -----------------------------
@dataclass
class ProjectPaths:
    root: Path
    pytools: Path
    apps: Path
    sim_data: Path
    sim_res: Path
    run_sie: Path
    vis: Path

    @staticmethod
    def discover(start: Path) -> "ProjectPaths":
        root = start.resolve()
        pytools = root / "pytools"
        apps = root / "apps"
        sim_data = root / "sim_data"
        sim_res = root / "sim_res"
        run_sie = root / "run_sie.py"
        vis = pytools / "visualization.py"
        return ProjectPaths(root, pytools, apps, sim_data, sim_res, run_sie, vis)

# ----------------------------- queued ProcRunner -----------------------------
class ProcRunner(QtCore.QObject):
    started  = QtCore.Signal(str)
    line     = QtCore.Signal(str)
    finished = QtCore.Signal(int)

    def __init__(self, cwd: Path | None = None, env: dict | None = None):
        super().__init__()
        self.cwd = str(cwd) if cwd else None
        self.env = os.environ.copy()
        if env:
            self.env.update(env)
        self._worker: threading.Thread | None = None
        self._stop_flag = threading.Event()
        self._proc: subprocess.Popen | None = None

    def run(self, args: list[str] | str):
        """Run a single job; reject if one is already running."""
        if isinstance(args, str):
            cmd = args
        else:
            cmd = " ".join(shlex.quote(x) for x in args)

        if self._worker and self._worker.is_alive():
            self.line.emit("[warn] A job is already running. Please cancel it first.\n")
            return

        self._stop_flag.clear()
        self._worker = threading.Thread(target=self._run_proc, args=(cmd,), daemon=True)
        self._worker.start()
        self.started.emit(cmd)

    def cancel_all(self, force_after: float = 3.0):
        """Terminate the running job and all its children (process group)."""
        self._stop_flag.set()
        if not self._proc or self._proc.poll() is not None:
            return
        try:
            # send SIGTERM to the whole group (works on POSIX)
            os.killpg(self._proc.pid, signal.SIGTERM)
        except Exception:
            try:
                self._proc.terminate()
            except Exception:
                pass
        # grace period
        t0 = time.time()
        while time.time() - t0 < force_after and self._proc.poll() is None:
            time.sleep(0.1)
        # if still alive, SIGKILL group
        if self._proc.poll() is None:
            try:
                os.killpg(self._proc.pid, signal.SIGKILL)
            except Exception:
                try:
                    self._proc.kill()
                except Exception:
                    pass

    # ---------------- internal ----------------
    def _run_proc(self, cmd: str):
        try:
            # Start the child in its own process group so we can kill the group.
            self._proc = subprocess.Popen(
                cmd,
                cwd=self.cwd,
                env=self.env,
                shell=True,
                bufsize=1,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                universal_newlines=True,
                preexec_fn=os.setsid  # new process group (POSIX)
            )
            assert self._proc.stdout is not None
            while True:
                ln = self._proc.stdout.readline()
                if ln:
                    self.line.emit(ln)
                else:
                    if self._proc.poll() is not None:
                        break
                if self._stop_flag.is_set():
                    # cancel_all() will do the killing; just break reading loop
                    break
            rc = self._proc.wait()
            self.finished.emit(rc)
        except FileNotFoundError as e:
            self.line.emit(f"[error] {e}\n")
            self.finished.emit(127)
        except Exception as e:
            self.line.emit(f"[error] {e}\n")
            self.finished.emit(1)
        finally:
            self._proc = None

# ----------------------------- Matplotlib widgets -----------------------------
class MplPane(QtWidgets.QWidget):
    def __init__(self, parent=None, with3d=False):
        super().__init__(parent)
        self.fig = Figure(figsize=(5, 4), tight_layout=True)
        if with3d:
            self.ax = self.fig.add_subplot(111, projection='3d')
        else:
            self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvas(self.fig)
        lay = QtWidgets.QVBoxLayout(self); lay.setContentsMargins(0,0,0,0)
        lay.addWidget(self.canvas)

    def clear(self):
        self.ax.cla()

    def draw(self):
        self.canvas.draw_idle()

# ----------------------------- main GUI -----------------------------
class HeliosGUI(QtWidgets.QMainWindow):
    def __init__(self, proj: ProjectPaths):
        super().__init__()
        self.setWindowTitle("HELIOS GUI")
        self.resize(1200, 800)
        self.proj = proj
        self.runner = ProcRunner(cwd=self.proj.root)
        self.runner.started.connect(self.on_started)
        self.runner.line.connect(self.on_line)
        self.runner.finished.connect(self.on_finished)

        self._make_ui()
        self._refresh_sim_list()

    # ------------------------- UI skeleton -------------------------
    def _make_ui(self):
        cw = QtWidgets.QWidget(); self.setCentralWidget(cw)
        main = QtWidgets.QVBoxLayout(cw)

        # top bar: project root + sim selector
        top = QtWidgets.QHBoxLayout()
        self.root_edit = QtWidgets.QLineEdit(str(self.proj.root))
        br_root = QtWidgets.QPushButton("Browse…")
        br_root.clicked.connect(self.pick_root)
        self.sim_combo = QtWidgets.QComboBox()
        self.refresh_btn = QtWidgets.QPushButton("↻")
        self.refresh_btn.setToolTip("Refresh simulations list")
        self.refresh_btn.clicked.connect(self._refresh_sim_list)
        self.new_sim_btn = QtWidgets.QPushButton("New sim…")
        self.new_sim_btn.clicked.connect(self.create_sim)

        # Open sim folder button with menu
        self.open_sim_btn = QtWidgets.QToolButton()
        self.open_sim_btn.setText("Open…")
        self.open_sim_btn.setPopupMode(QtWidgets.QToolButton.MenuButtonPopup)

        m = QtWidgets.QMenu(self.open_sim_btn)
        act_res = m.addAction("Open sim_res/<sim>")
        act_data = m.addAction("Open sim_data/<sim>")
        self.open_sim_btn.setMenu(m)

        def _open_res():
            sim = self.cur_sim()
            if not sim:
                QtWidgets.QMessageBox.warning(self, "No simulation", "Pick or create a simulation first.")
                return
            self._open_path(self.proj.sim_res / sim)

        def _open_data():
            sim = self.cur_sim()
            if not sim:
                QtWidgets.QMessageBox.warning(self, "No simulation", "Pick or create a simulation first.")
                return
            self._open_path(self.proj.sim_data / sim)

        act_res.triggered.connect(_open_res)
        act_data.triggered.connect(_open_data)

        logo_label = QtWidgets.QLabel()
        pm = QtGui.QPixmap(str(self.proj.root / "logo.png")).scaledToHeight(28, QtCore.Qt.SmoothTransformation)
        logo_label.setPixmap(pm)
        top.addWidget(logo_label, 0)

        # default click -> sim_res
        self.open_sim_btn.clicked.connect(_open_res)

        top.addWidget(self.open_sim_btn)
        top.addWidget(QtWidgets.QLabel("Project root:"))
        top.addWidget(self.root_edit, 1)
        top.addWidget(br_root)
        top.addSpacing(16)
        top.addWidget(QtWidgets.QLabel("Simulation:"))
        top.addWidget(self.sim_combo, 0)
        top.addWidget(self.refresh_btn)
        top.addWidget(self.new_sim_btn)
        main.addLayout(top)

        # Tabs
        tabs = QtWidgets.QTabWidget()
        tabs.addTab(self._tab_prepare(), "Pre-processing")
        tabs.addTab(self._tab_mesh(), "Mesh Viewer")
        tabs.addTab(self._tab_solve(), "Solver")
        tabs.addTab(self._tab_points(), "Generate Points")
        tabs.addTab(self._tab_post(), "Post-processing")
        tabs.addTab(self._tab_distributed(), "Distributed Run")
        tabs.addTab(self._tab_visualize(), "Visualization")
        tabs.addTab(self._tab_table(), "Generate Table")
        tabs.addTab(self._tab_settings(), "Settings")
        main.addWidget(tabs, 1)

        # Progress bar (global for long tasks)
        self.pb = QtWidgets.QProgressBar()
        self.pb.setTextVisible(True)
        self.pb.setFormat("Progress Bar")
        self.pb.setMaximum(1)   # 0/0 = busy, but Qt uses 0/0 for range; we force determinate later
        self.pb.setValue(0)
        main.addWidget(self.pb)

        # Log panel
        self.log = QtWidgets.QPlainTextEdit(); self.log.setReadOnly(True)
        self.log.setLineWrapMode(QtWidgets.QPlainTextEdit.NoWrap)
        self.log.setMaximumBlockCount(10000)
        main.addWidget(self.log, 1)

        # bottom row: Cancel + Clean output
        btn_layout = QtWidgets.QHBoxLayout()
        btn_layout.addStretch(1)

        # Clean output button
        self.clear_btn = QtWidgets.QPushButton("Clean console")
        self.clear_btn.setToolTip("Erase all text from the console above.")
        self.clear_btn.clicked.connect(self.clean_output)
        btn_layout.addWidget(self.clear_btn)

        # Cancel all button
        self.cancel_btn = QtWidgets.QPushButton("Cancel all")
        self.cancel_btn.setToolTip("Stop all running jobs.")
        self.cancel_btn.clicked.connect(self.runner.cancel_all)
        btn_layout.addWidget(self.cancel_btn)

        main.addLayout(btn_layout)

        # Status
        self.status = self.statusBar()

    # ------------------------- helper widgets -------------------------
    def _num(self, label, val, minv=None, maxv=None, step=1):
        w = QtWidgets.QSpinBox()
        if minv is not None: w.setMinimum(minv)
        if maxv is not None: w.setMaximum(maxv)
        w.setValue(val); w.setSingleStep(step)
        box = QtWidgets.QWidget(); lay = QtWidgets.QHBoxLayout(box); lay.setContentsMargins(0,0,0,0)
        lay.addWidget(QtWidgets.QLabel(label)); lay.addWidget(w)
        return box, w

    def _text(self, label, val=""):
        e = QtWidgets.QLineEdit(val)
        box = QtWidgets.QWidget(); lay = QtWidgets.QHBoxLayout(box); lay.setContentsMargins(0,0,0,0)
        lay.addWidget(QtWidgets.QLabel(label)); lay.addWidget(e)
        return box, e

    def _check(self, label, val=False):
        c = QtWidgets.QCheckBox(label); c.setChecked(val)
        return c

    def _file_row(self, label, mode="file"):
        e = QtWidgets.QLineEdit()
        b = QtWidgets.QPushButton("…")
        def pick():
            if mode == "file":
                p, _ = QtWidgets.QFileDialog.getOpenFileName(self, label, str(self.proj.root))
            else:
                p = QtWidgets.QFileDialog.getExistingDirectory(self, label, str(self.proj.root))
            if p:
                e.setText(p)
        b.clicked.connect(pick)
        box = QtWidgets.QWidget(); lay = QtWidgets.QHBoxLayout(box); lay.setContentsMargins(0,0,0,0)
        lay.addWidget(QtWidgets.QLabel(label)); lay.addWidget(e, 1); lay.addWidget(b)
        return box, e
    
    def clean_output(self):
        """Clean console and reset progress bar to Ready."""
        self.log.clear()
        self.status.showMessage("Console cleared.")
        self.pb.setRange(0, 1)
        self.pb.setValue(1)
        self.pb.setFormat("Ready")

    # ---------- progress context (sum-of-totals, count-of-dones) ----------
    def _pg_reset(self, label: str, tag: str):
        """
        tag: 'solve' or 'post'. We ignore any other tags (e.g., 'dist').
        """
        self._pg = {"label": label, "tag": tag, "total": 0, "done": 0}
        # start busy until we see at least one 'total'
        self.pb.setRange(0, 0)
        self.pb.setFormat(label)

    def _pg_set_range(self):
        # switch to determinate once we have a positive total
        if self._pg["total"] > 0:
            self.pb.setRange(0, self._pg["total"])
            # clamp done to total
            v = min(self._pg["done"], self._pg["total"])
            self.pb.setValue(v)
            self.pb.setFormat(f'{self._pg["label"]} %p%')

    def _progress_parse_line(self, s: str):
        """
        Implements:
        Total progress = sum of all '[<tag>-progress] total N' (tag in {'solve','post'})
        Ongoing progress = count of appearances of '[<tag>-progress] done K'
        Ignore any 'dist-progress' lines entirely.
        """
        tag = self._pg.get("tag") if hasattr(self, "_pg") else None
        if tag not in ("solve", "post"):
            return  # we're not in a solve/post context; ignore

        # Match "[solve-progress] total N" or "[post-progress] total N"
        m = re.search(rf"\[{tag}-progress\]\s+total\s+(\d+)", s)
        if m:
            self._pg["total"] += int(m.group(1))
            self._pg_set_range()  # may flip from busy to determinate
            return

        # Match "[solve-progress] done K" or "[post-progress] done K"
        if re.search(rf"\[{tag}-progress\]\s+done\s+\d+", s):
            self._pg["done"] += 1  # each appearance counts as +1
            if self.pb.maximum() > 0:  # only if we already know total
                v = min(self._pg["done"], self._pg["total"])
                self.pb.setValue(v)
            return

    # ------------------------- tabs -------------------------
    def _tab_prepare(self):
        w = QtWidgets.QWidget(); f = QtWidgets.QFormLayout(w)
        self.mode = QtWidgets.QComboBox(); 
        self.mode.addItems(["isolated scatterer(s) in homogeneous background",
                            "periodic structure in homogeneous background",
                            "isolated scatterer(s) in layered background"]); 
        self.mode.setCurrentIndex(2)
        # mesh source (optional): .mphtxt or Gmsh .msh (v4 ASCII recommended)
        self.mesh_row, self.mesh_edit = self._file_row("Mesh file (.mphtxt / .msh):", mode="file")
        self.overwrite = self._check("Overwrite existing folder")
        self.mat_cols = QtWidgets.QComboBox()
        self.mat_cols.addItems(["ε (Re, Im)", "n,k (convert to ε)"])
        self.mat_cols.setToolTip("How to interpret materials tables for jobwriter via run_sie.py --materials")
        self.spline_interp = self._check("Use spline (not linear) interpolation for materials", True)
        f.addRow("Mode:", self.mode)
        f.addRow(self.mesh_row)
        f.addRow(self.overwrite)
        f.addRow("Materials columns:", self.mat_cols)
        f.addRow(self.spline_interp)
        self.btn_prep = QtWidgets.QPushButton("Prepare simulation")
        self.btn_prep.clicked.connect(self.do_prepare)
        f.addRow(self.btn_prep)
        return w

    def _tab_solve(self):
        w = QtWidgets.QWidget(); g = QtWidgets.QGridLayout(w)
        self.acc = self._check("Enable increased accuracy")
        level_row, self.level = self._num("Linear system solver [1 = iterative (CG), 2 = direct (LU)]:", 2, 0, 9)
        th_row, self.threads = self._num("Threads (0 = max available):", 0, 0, 256)
        self.jobs_row, self.jobs = self._text("Jobs' selection:")
        self.jobs.setPlaceholderText("e.g. 1,3,5-8")
        self.sel_by_lambda = self._check("Select by λ [nm] instead of jobs")
        self.lambda_row, self.lambda_edit = self._text("λ list [nm]:")
        self.lambda_edit.setPlaceholderText("e.g. 500, 650, 750")    
        self._bind_lambda_toggle(self.sel_by_lambda, self.jobs_row, self.lambda_row)
        etm_row, self.etm = self._num("Evaluation terms (only periodic, for isolated put 0):", 0, 0, 64)
        self.btn_solve = QtWidgets.QPushButton("Run solver")
        self.btn_solve.clicked.connect(self.do_solve)
        row=0
        for w2 in [self.acc, level_row, th_row, self.sel_by_lambda, self.jobs_row, self.lambda_row, etm_row]:
            g.addWidget(w2, row, 0, 1, 2); row+=1
        g.addWidget(self.btn_solve, row, 0, 1, 2)
        return w

    def _tab_post(self):
        w = QtWidgets.QWidget(); g = QtWidgets.QGridLayout(w)
        self.acc_post = self._check("Enable increased accuracy")
        thp_row, self.th_post = self._num("Threads (0 = max available):", 0, 0, 256)
        self.jobs_post_row, self.jobs_post = self._text("Jobs' selection:")
        self.jobs_post.setPlaceholderText("e.g. 1,3,5-8")
        self.sel_by_lambda_post = self._check("Select by λ [nm] instead of jobs")
        self.lambda_post_row, self.lambda_post_edit = self._text("λ list [nm]:")
        self.lambda_post_edit.setPlaceholderText("e.g. 500, 650, 750")
        self._bind_lambda_toggle(self.sel_by_lambda_post, self.jobs_post_row, self.lambda_post_row)
        etm_row, self.etm_post = self._num("Evaluation terms (only periodic, for isolated put 0):", 0, 0, 64)
        self.points_row, self.points_edit = self._text("Points file name:")
        self.points_edit.setPlaceholderText("e.g. points_file or points_file.pos")
        self.btn_post = QtWidgets.QPushButton("Run post-processing")
        self.btn_post.clicked.connect(self.do_post)
        row=0
        for w2 in [self.acc_post, thp_row, etm_row, self.sel_by_lambda_post, 
                   self.jobs_post_row, self.lambda_post_row, self.points_row]:
            g.addWidget(w2, row, 0, 1, 2); row+=1
        g.addWidget(self.btn_post, row, 0, 1, 2)
        return w

    def _tab_points(self):
        def L(html: str) -> QtWidgets.QLabel:
            lab = QtWidgets.QLabel(html)
            lab.setTextFormat(Qt.RichText)  # allow <sub>
            return lab

        w = QtWidgets.QWidget()
        form = QtWidgets.QFormLayout(w)
        form.setLabelAlignment(Qt.AlignRight)
        form.setFormAlignment(Qt.AlignLeft)
        form.setHorizontalSpacing(10)
        form.setVerticalSpacing(4)
        form.setContentsMargins(10, 6, 10, 6)

        # Plane selector
        self.plane = QtWidgets.QComboBox()
        self.plane.addItems(["xz", "xy", "yz"])
        form.addRow(L("Plane:"), self.plane)

        # Create edits directly (no container helper)
        self.xmin = QtWidgets.QLineEdit(); self.xmax = QtWidgets.QLineEdit()
        self.ymin = QtWidgets.QLineEdit(); self.ymax = QtWidgets.QLineEdit()
        self.zmin = QtWidgets.QLineEdit(); self.zmax = QtWidgets.QLineEdit()

        self.step  = QtWidgets.QLineEdit()
        self.stepx = QtWidgets.QLineEdit()
        self.stepy = QtWidgets.QLineEdit()
        self.stepz = QtWidgets.QLineEdit()

        self.x0 = QtWidgets.QLineEdit()
        self.y0 = QtWidgets.QLineEdit()
        self.z0 = QtWidgets.QLineEdit()

        # examples for multi-values in GUI
        self.x0.setPlaceholderText("e.g. -200, 0, 200")
        self.y0.setPlaceholderText("e.g. 0, 100")
        self.z0.setPlaceholderText("e.g. -500, 500")

        self.outname = QtWidgets.QLineEdit()
        self.outname.setPlaceholderText("e.g. points_file1 or points_file1.pos")

        # Keep label widgets so we can hide/show both label and editor
        self._row_widgets = {
            "xmin": (L("x<sub>min</sub>:"), self.xmin),
            "xmax": (L("x<sub>max</sub>:"), self.xmax),
            "ymin": (L("y<sub>min</sub>:"), self.ymin),
            "ymax": (L("y<sub>max</sub>:"), self.ymax),
            "zmin": (L("z<sub>min</sub>:"), self.zmin),
            "zmax": (L("z<sub>max</sub>:"), self.zmax),

            "step":  (L("step (uniform):"), self.step),
            "stepx": (L("step<sub>x</sub>:"), self.stepx),
            "stepy": (L("step<sub>y</sub>:"), self.stepy),
            "stepz": (L("step<sub>z</sub>:"), self.stepz),

            "x0": (L("x<sub>0</sub>:"), self.x0),
            "y0": (L("y<sub>0</sub>:"), self.y0),
            "z0": (L("z<sub>0</sub>:"), self.z0),

            "outname": (L("Points file name:"), self.outname),
        }

        # Add all rows once; we will toggle visibility later
        for key in ["xmin","xmax","ymin","ymax","zmin","zmax",
                    "step","stepx","stepy","stepz",
                    "x0","y0","z0","outname"]:
            lab, edt = self._row_widgets[key]
            form.addRow(lab, edt)

        # Generate button
        self.btn_pts = QtWidgets.QPushButton("Generate points file")
        self.btn_pts.clicked.connect(self.do_points)
        form.addRow(self.btn_pts)

        # Wire behaviors
        self.plane.currentTextChanged.connect(self._points_update_visibility)
        for e in (self.stepx, self.stepy, self.stepz):
            e.textChanged.connect(self._points_update_steps)
        self.step.textChanged.connect(self._points_update_steps)

        # Initial state
        self._points_update_visibility()
        self._points_update_steps()
        return w

    def _tab_table(self):
        w = QtWidgets.QWidget()
        f = QtWidgets.QFormLayout(w)

        # Filename and text/binary mode
        self.tab_filename_row, self.tab_filename = self._text("Output filename:")
        self.tab_filename.setText("erfc.bin")
        self.tab_text_mode = self._check("Save as text (instead of binary)", False)

        # Numeric inputs for table bounds and increments
        self.tab_maxRe = QtWidgets.QDoubleSpinBox()
        self.tab_maxRe.setRange(0.1, 1000.0); self.tab_maxRe.setValue(10.0); self.tab_maxRe.setSingleStep(1.0)
        
        self.tab_maxIm = QtWidgets.QDoubleSpinBox()
        self.tab_maxIm.setRange(0.1, 1000.0); self.tab_maxIm.setValue(10.0); self.tab_maxIm.setSingleStep(1.0)
        
        self.tab_incRe = QtWidgets.QDoubleSpinBox()
        self.tab_incRe.setRange(0.001, 1.0); self.tab_incRe.setValue(0.1); self.tab_incRe.setSingleStep(0.01); self.tab_incRe.setDecimals(3)
        
        self.tab_incIm = QtWidgets.QDoubleSpinBox()
        self.tab_incIm.setRange(0.001, 1.0); self.tab_incIm.setValue(0.1); self.tab_incIm.setSingleStep(0.01); self.tab_incIm.setDecimals(3)

        f.addRow(self.tab_filename_row)
        f.addRow(self.tab_text_mode)
        f.addRow("Max Re:", self.tab_maxRe)
        f.addRow("Max Im:", self.tab_maxIm)
        f.addRow("Inc Re:", self.tab_incRe)
        f.addRow("Inc Im:", self.tab_incIm)

        # Action button
        self.btn_table = QtWidgets.QPushButton("Generate Table")
        self.btn_table.clicked.connect(self.do_table)
        f.addRow(self.btn_table)

        return w

    def _points_update_visibility(self):
        """Show only fields needed for the chosen plane."""
        plane = self.plane.currentText()

        # Start by hiding all
        def hide(keys):
            for k in keys:
                lab, edt = self._row_widgets[k]
                lab.setVisible(False); edt.setVisible(False)

        def show(keys):
            for k in keys:
                lab, edt = self._row_widgets[k]
                lab.setVisible(True); edt.setVisible(True)

        all_keys = ["xmin","xmax","ymin","ymax","zmin","zmax","x0","y0","z0"]
        hide(all_keys)

        if plane == "xz":
            show(["xmin","xmax","zmin","zmax","y0"])
        elif plane == "xy":
            show(["xmin","xmax","ymin","ymax","z0"])
        else:  # yz
            show(["ymin","ymax","zmin","zmax","x0"])

        # Steps + output name are always visible
        show(["step","stepx","stepy","stepz","outname"])

    def _points_update_steps(self):
        """If uniform step is set, lock axis steps; if any axis step is set, lock uniform."""
        has_uniform = self.step.text().strip() != ""
        has_axis = any(e.text().strip() != "" for e in (self.stepx, self.stepy, self.stepz))

        # Uniform typed -> disable axis edits
        for e in (self.stepx, self.stepy, self.stepz):
            e.setEnabled(not has_uniform)

        # Any axis typed -> disable uniform edit
        self.step.setEnabled(not has_axis)

    def resizeEvent(self, ev):
        super().resizeEvent(ev)
        try:
            # keep the Options panel at a fixed fraction of current window height
            h = int(self.height() * getattr(self, "viz_height_ratio", 0.55))
            if hasattr(self, "viz_ctrl") and self.viz_ctrl is not None:
                self.viz_ctrl.setMinimumHeight(h)
                self.viz_ctrl.setMaximumHeight(h)
        except Exception:
            pass        

    def _tab_visualize(self):
        """
        Visualization tab with three modes: Cross sections, Mie solution, Fields.
        Only the controls relevant to the chosen mode are shown.
        """
        w = QtWidgets.QWidget()
        outer = QtWidgets.QHBoxLayout(w)

        # ---- LEFT: controls panel ----
        ctrl = QtWidgets.QGroupBox("Options")
        form = QtWidgets.QFormLayout(ctrl)
        self.viz_ctrl = ctrl  # remember for height control
        form.setContentsMargins(6, 6, 6, 6)
        form.setHorizontalSpacing(6)
        form.setVerticalSpacing(6)

        # --- Top row: visible mode selectors (mutually exclusive) ---
        self.v_sel_xs  = QtWidgets.QRadioButton("Cross sections")
        self.v_sel_nf  = QtWidgets.QRadioButton("Fields")
        self.v_sel_rt  = QtWidgets.QRadioButton("R/T")
        self.v_sel_nf.setChecked(True)  # default: Fields

        sel_grp = QtWidgets.QButtonGroup(self)
        sel_grp.setExclusive(True)
        sel_grp.addButton(self.v_sel_xs)
        sel_grp.addButton(self.v_sel_nf)
        sel_grp.addButton(self.v_sel_rt)

        top_row = QtWidgets.QHBoxLayout()
        top_row.setContentsMargins(0, 0, 0, 0)
        top_row.setSpacing(8)
        top_row.addWidget(self.v_sel_xs)
        top_row.addWidget(self.v_sel_nf)
        top_row.addWidget(self.v_sel_rt)
        top_row.addStretch(1)
        form.addRow("Plot type:", top_row)

        # --- Cross-section series (QListWidget) ---
        self.v_series = QtWidgets.QListWidget()
        self.v_series.addItems([
            "scattering cross  section (scs)",
            "absorption cross  section (acs)",
            "extinction cross  section (ext)",
        ])
        self.v_series.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)
        for i in range(self.v_series.count()):
            self.v_series.item(i).setSelected(True)

        series_row = QtWidgets.QWidget()
        series_lay = QtWidgets.QHBoxLayout(series_row)
        series_lay.setContentsMargins(0, 0, 0, 0)
        series_lay.addWidget(self.v_series)
        form.addRow(series_row)

        # --- Simulations selector (multi-select) ---
        self.v_simlist = QtWidgets.QListWidget()
        self.v_simlist.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAsNeeded)
        self.v_simlist.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.v_simlist.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)
        self.v_simlist.setMinimumHeight(80)
        self.v_simlist.setMaximumHeight(80)
        self.v_simlist.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)
        self.v_simlist.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.v_simlist.setFocusPolicy(QtCore.Qt.WheelFocus)
        self.v_simlist.setVerticalScrollMode(QtWidgets.QAbstractItemView.ScrollPerPixel)
        sim_row = QtWidgets.QWidget()
        sim_lay = QtWidgets.QVBoxLayout(sim_row)
        sim_lay.setContentsMargins(0, 0, 0, 0)
        sim_lay.addWidget(self.v_simlist)
        self.v_sim_label = QtWidgets.QLabel("Simulations:")
        form.addRow(self.v_sim_label, sim_row)

        # --- Mie parameters ---
        mie_row = QtWidgets.QWidget()
        mie_lay = QtWidgets.QGridLayout(mie_row)
        mie_lay.setContentsMargins(0, 0, 0, 0)
        mie_lay.setHorizontalSpacing(6)
        mie_lay.setVerticalSpacing(3)

        self.v_mie_R = QtWidgets.QDoubleSpinBox(); self.v_mie_R.setRange(0.0, 1e9); self.v_mie_R.setDecimals(3); self.v_mie_R.setValue(50.0)
        self.v_mie_mat = QtWidgets.QLineEdit(); self.v_mie_mat.setPlaceholderText("eps=Re,Im  |  nk=n,k  |  @file  |  Au")
        self.v_mie_bg  = QtWidgets.QLineEdit(); self.v_mie_bg.setPlaceholderText("eps=1,0  (vacuum default)")

        mie_lay.addWidget(QtWidgets.QLabel("R [nm]:"),   0, 0); mie_lay.addWidget(self.v_mie_R,  0, 1)
        mie_lay.addWidget(QtWidgets.QLabel("Material:"), 1, 0); mie_lay.addWidget(self.v_mie_mat,1, 1)
        mie_lay.addWidget(QtWidgets.QLabel("Background:"),2, 0); mie_lay.addWidget(self.v_mie_bg, 2, 1)
        self.v_mie_enable = QtWidgets.QPushButton("Show")
        self.v_mie_enable.toggled.connect(lambda on: self.v_mie_enable.setText("Hide" if on else "Show"))        
        self.v_mie_enable.setCheckable(True)
        mie_hdr = QtWidgets.QWidget()
        mie_hdr_lay = QtWidgets.QFormLayout(mie_hdr)
        mie_hdr_lay.setContentsMargins(0, 0, 0, 0)
        mie_hdr_lay.setHorizontalSpacing(6)
        mie_hdr_lay.setVerticalSpacing(3)
        mie_hdr_lay.addRow("Mie solution:", self.v_mie_enable)
        form.addRow(mie_hdr)
        form.addRow(mie_row)

        # --- Near-field controls (only for "Fields") ---
        nf_row = QtWidgets.QWidget()
        nfl = QtWidgets.QFormLayout(nf_row)
        nfl.setContentsMargins(0, 0, 0, 0)
        nfl.setHorizontalSpacing(6)
        nfl.setVerticalSpacing(3)

        self.v_points = QtWidgets.QLineEdit(); self.v_points.setPlaceholderText("filename.pos (e.g. points_file1.pos)")
        self.v_field  = QtWidgets.QComboBox(); self.v_field.addItems(["e","h"])
        self.v_part   = QtWidgets.QComboBox(); self.v_part.addItems(["inc","sc","tot"])
        self.v_func   = QtWidgets.QComboBox(); self.v_func.addItems(["intensity", "magnitude", "real", "imag", "phase"])
        self.v_quan   = QtWidgets.QComboBox(); self.v_quan.addItems(["total", "x", "y", "z"])
        self.v_mode   = QtWidgets.QComboBox(); self.v_mode.addItems(["2d","3d","1d"])
        self.v_line   = QtWidgets.QLineEdit(); self.v_line.setPlaceholderText("x=0 or y=... or z=...")
        self.v_lambda = QtWidgets.QDoubleSpinBox(); self.v_lambda.setDecimals(3); self.v_lambda.setRange(0.0, 1e9); self.v_lambda.setValue(0.0)
        self.v_lambda.setToolTip("If > 0, choose closest wavelength [nm]")
        self.v_subjob = QtWidgets.QSpinBox()
        self.v_subjob.setRange(1, 999999)
        self.v_subjob.setValue(1)
        self.v_show_ints = QtWidgets.QPushButton("Show")
        self.v_show_ints.toggled.connect(lambda on: self.v_show_ints.setText("Hide" if on else "Show"))
        self.v_show_ints.setCheckable(True)
        self.v_log = QtWidgets.QPushButton("Enable")
        self.v_log.toggled.connect(lambda on: self.v_log.setText("Disable" if on else "Enable"))
        self.v_log.setCheckable(True)

        nfl.addRow("Points name:", self.v_points)
        nfl.addRow("Field type:", self.v_field)
        nfl.addRow("Field part:", self.v_part)
        nfl.addRow("Plot function:", self.v_func)
        nfl.addRow("Plot quantity:", self.v_quan)
        nfl.addRow("Plot dimension:", self.v_mode)
        nfl.addRow("1D line:", self.v_line)
        nfl.addRow("λ [nm]:", self.v_lambda)
        nfl.addRow("Interfaces:", self.v_show_ints)
        nfl.addRow("Logscale:", self.v_log)

        form.addRow(nf_row)

        # R/T controls
        rt_row = QtWidgets.QWidget()
        rtl = QtWidgets.QFormLayout(rt_row)
        rtl.setContentsMargins(0, 0, 0, 0)
        rtl.setHorizontalSpacing(6)
        rtl.setVerticalSpacing(3)

        self.v_rt_points = QtWidgets.QLineEdit(); self.v_rt_points.setPlaceholderText("filename.pos (e.g. points_file1.pos)")
        self.v_rt_prop = QtWidgets.QComboBox(); self.v_rt_prop.addItems(["auto","+z","-z"])   

        rtl.addRow("Points:", self.v_rt_points)
        rtl.addRow("Propagation:", self.v_rt_prop)

        form.addRow(rt_row)

        # --- Subjob ---
        form.addRow("Subjob:", self.v_subjob)

        # --- Save / file output ---
        self.v_save = self._check("Save to file")
        self.v_saverow, self.v_savepath = self._file_row("Filename: ", mode="file")
        self.v_saverow.setEnabled(False)
        self.v_save.toggled.connect(lambda s: self.v_saverow.setEnabled(s))
        form.addRow(self.v_save)
        form.addRow(self.v_saverow)

        # --- Render button ---
        self.btn_vis = QtWidgets.QPushButton("Render")
        self.btn_vis.clicked.connect(self.do_visualize)
        form.addRow(self.btn_vis)

        # --- Visibility logic ---
        def _apply_mode():
            is_xs  = self.v_sel_xs.isChecked()
            is_nf  = self.v_sel_nf.isChecked()
            is_rt  = self.v_sel_rt.isChecked()

            series_row.setVisible(is_xs)
            sim_row.setVisible(is_xs)
            self.v_sim_label.setVisible(is_xs)
            mie_hdr.setVisible(is_xs)
            mie_row.setVisible(is_xs)
            nf_row.setVisible(is_nf)
            rt_row.setVisible(is_rt)

        self.v_sel_xs.toggled.connect(_apply_mode)
        self.v_sel_nf.toggled.connect(_apply_mode)
        self.v_sel_rt.toggled.connect(_apply_mode)
        _apply_mode()  # initialize

        outer.addWidget(ctrl, 0)

        # ---- RIGHT: preview canvas ----
        self.v_plot = MplPane(with3d=False)
        outer.addWidget(self.v_plot, 1)
        return w

    def _tab_mesh(self):
        w = QtWidgets.QWidget()
        outer = QtWidgets.QHBoxLayout(w)

        # -------- Left options --------
        left = QtWidgets.QGroupBox("Options")
        form = QtWidgets.QFormLayout(left)

        # mode chooser
        self.m_mode = QtWidgets.QComboBox()
        self.m_mode.addItems(["3D", "2D: XY", "2D: XZ", "2D: YZ", "1D: z-interfaces"])

        # your existing widgets
        self.m_decimate = QtWidgets.QSpinBox()
        self.m_decimate.setRange(1, 1000)
        self.m_decimate.setValue(1)
        self.m_alpha = QtWidgets.QDoubleSpinBox()
        self.m_alpha.setRange(0.05,1.0)
        self.m_alpha.setSingleStep(0.05)
        self.m_alpha.setValue(1.00)

        # Save to file
        self.m_save = self._check("Save to file")
        self.m_saverow, self.m_savepath = self._file_row("Filename:", mode="file")
        self.m_saverow.setEnabled(False)
        self.m_save.toggled.connect(lambda s: self.m_saverow.setEnabled(s))

        self.btn_mesh = QtWidgets.QPushButton("Load / Render")
        self.btn_mesh.clicked.connect(self.do_mesh)

        form.addRow("Mode:", self.m_mode)
        form.addRow("Decimate:", self.m_decimate)
        form.addRow("Alpha:", self.m_alpha)
        form.addRow(self.m_save)
        form.addRow(self.m_saverow)
        form.addRow(self.btn_mesh)
        outer.addWidget(left, 0)

        # -------- Right: Matplotlib canvas (2D or 3D) --------
        self.m_fig = Figure(figsize=(9, 7))
        # We’ll swap projection depending on mode; start with 3D:
        self.m_ax = self.m_fig.add_subplot(111, projection="3d")
        self.m_canvas = FigureCanvas(self.m_fig)

        right = QtWidgets.QVBoxLayout()
        right.addWidget(self.m_canvas)
        outer.addLayout(right, 1)

        return w

    def _tab_distributed(self):
        w = QtWidgets.QWidget(); g = QtWidgets.QGridLayout(w)
        # Which launcher script
        self.dist_mode = QtWidgets.QComboBox()
        self.dist_mode.addItems(["isolated scatterer(s) in homogeneous background",
                                 "periodic structure in homogeneous background",
                                 "isolated scatterer(s) in layered background"])
        self.dist_mode.setCurrentIndex(2)                                
        # Nodes: names or count
        self.dist_nodes_row, self.dist_nodes = self._text("Nodes:")
        self.dist_nodes.setPlaceholderText("e.g. n01,n03,n04")
        # Subcommand & options
        self.dist_subcmd = QtWidgets.QComboBox(); self.dist_subcmd.addItems(["Run solver","Run post-processing"])
        lvl_row, self.dist_level = self._num("Linear system solver [1 = iterative (CG), 2 = direct (LU)]:", 2, 0, 9)
        th_row, self.dist_threads = self._num("Threads (0 = max available):", 0, 0, 256)
        etm_row, self.dist_etm = self._num("Evaluation terms (only periodic, for isolated put 0):", 0, 0, 64)
        self.dist_acc = self._check("Enable increased accuracy")
        self.dist_jobs_row, self.dist_jobs = self._text("Jobs' selection:")
        self.dist_jobs.setPlaceholderText("e.g. 1,3,5-8")
        self.dist_sel_by_lambda = self._check("Select by λ [nm] instead of jobs")
        self.dist_lambda_row, self.dist_lambda_edit = self._text("λ list [nm]:")
        self.dist_lambda_edit.setPlaceholderText("e.g. 500, 650, 750")
        self._bind_lambda_toggle(self.dist_sel_by_lambda, self.dist_jobs_row, self.dist_lambda_row)
        self.dist_points_row, self.dist_points = self._text("Points file name for post-processing:")
        self.dist_points.setPlaceholderText("e.g. points_file1 or points_file1.pos")
        # Run button
        self.btn_dist = QtWidgets.QPushButton("Run distributed simulation")
        self.btn_dist.clicked.connect(self.do_distributed)
        # Layout
        row=0
        for lab, wid in [("Simulation type", self.dist_mode), ("Simulation mode", self.dist_subcmd)]:
            g.addWidget(QtWidgets.QLabel(lab+":"), row, 0); g.addWidget(wid, row, 1); row+=1
        for rw in [self.dist_nodes_row, self.dist_acc, lvl_row, th_row, etm_row,
                   self.dist_sel_by_lambda, self.dist_jobs_row, self.dist_lambda_row,
                   self.dist_points_row, self.btn_dist]:
            g.addWidget(rw, row, 0, 1, 2); row += 1
        return w

    # ------------------------- actions -------------------------
    def pick_root(self):
        p = QtWidgets.QFileDialog.getExistingDirectory(self, "Pick project root", str(self.proj.root))
        if p:
            self.root_edit.setText(p)
            self.proj = ProjectPaths.discover(Path(p))
            self.run_sie_edit.setText(str(self.proj.run_sie))
            self.vis_edit.setText(str(self.proj.vis))
            self._refresh_sim_list()

    def _open_path(self, p: Path):
        """Open a directory in a Qt file dialog where you can view and open files."""
        if not p.exists():
            QtWidgets.QMessageBox.warning(self, "Not found", f"{p}\ndoes not exist.")
            return

        dlg = QtWidgets.QFileDialog(self, f"Browse: {p}")
        dlg.setFileMode(QtWidgets.QFileDialog.ExistingFile)  # show files & folders
        dlg.setOption(QtWidgets.QFileDialog.DontUseNativeDialog, True)
        dlg.setDirectory(str(p))
        dlg.setWindowModality(Qt.ApplicationModal)

        # Optional: allow double-click to open with the default system application
        def _open_selected():
            files = dlg.selectedFiles()
            if not files:
                return
            for f in files:
                QDesktopServices.openUrl(QUrl.fromLocalFile(f))

        dlg.accepted.connect(_open_selected)
        dlg.exec()

    def create_sim(self):
        name, ok = QtWidgets.QInputDialog.getText(self, "New simulation", "Name (folder under sim_data/):")
        if ok and name.strip():
            (self.proj.sim_data / name).mkdir(parents=True, exist_ok=True)
            (self.proj.sim_res / name).mkdir(parents=True, exist_ok=True)
            self._refresh_sim_list(select=name)

    def _refresh_sim_list(self, select: str | None = None):
        self.sim_combo.clear()
        sims = []
        if self.proj.sim_data.exists():
            sims = [p.name for p in sorted(self.proj.sim_data.iterdir()) if p.is_dir()]
        self.sim_combo.addItems(sims)
        # Also refresh Visualization tab multi-sim list (if it exists)
        if hasattr(self, "v_simlist") and self.v_simlist is not None:
            cur = self.cur_sim()
            self.v_simlist.clear()
            self.v_simlist.addItems(sims)
            # default selection: current sim (if present)
            for i in range(self.v_simlist.count()):
                it = self.v_simlist.item(i)
                if it.text() == cur:
                    it.setSelected(True)

    def cur_sim(self) -> str:
        return self.sim_combo.currentText().strip()

    def _collect_env(self) -> dict:
        env = {}
        for r in range(self.env_table.rowCount()):
            k = self.env_table.item(r,0).text() if self.env_table.item(r,0) else ""
            v = self.env_table.item(r,1).text() if self.env_table.item(r,1) else ""
            if k:
                env[k]=v
        return env

    def add_env_row(self):
        r = self.env_table.rowCount(); self.env_table.insertRow(r)
        self.env_table.setItem(r, 0, QtWidgets.QTableWidgetItem(""))
        self.env_table.setItem(r, 1, QtWidgets.QTableWidgetItem(""))

    def on_started(self, cmd: str):
        self.status.showMessage("Running: " + cmd)
        self.log.appendPlainText("$ " + cmd)
        if " post " in cmd:
            self._pg_reset("Post-processing… ", "post")
        elif " solve " in cmd:
            self._pg_reset("Solving… ", "solve")
        else:
            # not a solve/post: show busy but we won't parse progress
            self._pg_reset("Working… ", "none")

    def on_line(self, s: str):
        self.log.appendPlainText(s.rstrip("\n"))
        self._progress_parse_line(s)

    def on_finished(self, rc: int):
        self.status.showMessage(f"Done (rc={rc})")
        self.log.appendPlainText(f"[exit {rc}]\n")
        # Freeze bar to final state (Done/Failed), then move to 'Ready'
        self.pb.setRange(0, 1)
        self.pb.setValue(1 if rc == 0 else 0)
        self.pb.setFormat("Done" if rc == 0 else "Failed")
        # final ready state
        self.pb.setRange(0, 1)
        self.pb.setValue(1)
        self.pb.setFormat("Ready")

    # ------------------------- subprocess wrappers -------------------------
    def py(self) -> str:
        return shlex.quote(self.py_exec.text().strip() or sys.executable)

    def runsie(self) -> str:
        return shlex.quote(self.run_sie_edit.text().strip())

    def vispy(self) -> str:
        return shlex.quote(self.vis_edit.text().strip())

    def do_prepare(self):
        sim = self.cur_sim()
        if not sim:
            QtWidgets.QMessageBox.warning(self, "No simulation", "Pick or create a simulation first.")
            return

        mode_idx = self.mode.currentIndex()
        cmd = [self.py(), self.runsie(), "prepare", sim, "--mode", str(mode_idx)]
        mesh_path = self.mesh_edit.text().strip()
        if mesh_path: cmd += ["--meshfile", mesh_path]
        mat = "nk" if self.mat_cols.currentIndex() == 1 else "eps"
        cmd += ["--materials", mat]
        if self.spline_interp.isChecked(): cmd += ["--spline"]
        # if sim_res/<sim> exists and Overwrite is NOT checked, ask to overwrite
        out_dir = self.proj.sim_res / sim
        want_overwrite = self.overwrite.isChecked()
        if out_dir.exists() and not want_overwrite:
            btn = QtWidgets.QMessageBox.question(
                self, "Overwrite existing outputs?",
                f"Folder already exists:\n{out_dir}\n\nOverwrite it?",
                QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No,
                QtWidgets.QMessageBox.No
            )
            if btn == QtWidgets.QMessageBox.Yes:
                want_overwrite = True

        if want_overwrite:
            cmd += ["--overwrite"]

        self.runner.env.update(self._collect_env())
        self.runner.run(cmd)

    def do_solve(self):
        sim = self.cur_sim()
        if not sim:
            QtWidgets.QMessageBox.warning(self, "No simulation", "Pick or create a simulation first.")
            return
        mode_idx = self.mode.currentIndex()
        cmd = [self.py(), self.runsie(), "solve", sim,
               "--mode", str(mode_idx),
               "-l", str(self.level.value()), "-th", str(self.threads.value())]
        if self.acc.isChecked(): cmd += ["-a"]
        jobs_spec = self.jobs.text().strip()
        if self.sel_by_lambda.isChecked():
            lam_csv = self.lambda_edit.text().strip()
            if lam_csv:
                cmd += ["--lambdas", lam_csv]   # let run_sie.py do λ→jobs
        elif jobs_spec:
            cmd += ["-j", jobs_spec]
        if self.etm.value() != 0: cmd += ["--etm", str(self.etm.value())]
        self.runner.env.update(self._collect_env())
        self.runner.run(cmd)

    def do_post(self):
        sim = self.cur_sim()
        if not sim:
            QtWidgets.QMessageBox.warning(self, "No simulation", "Pick or create a simulation first.")
            return

        pts_raw = self.points_edit.text().strip().split()
        if not pts_raw:
            QtWidgets.QMessageBox.warning(self, "No points", "Provide one or more point-set names.")
            return

        # Resolve stems like 'points_xz' to full paths
        resolved = []
        base = self.proj.root / "sim_res" / sim / "points"
        for t in pts_raw:
            tpath = Path(t)
            if tpath.is_absolute() and tpath.exists():
                resolved.append(str(tpath))
                continue
            # Look under sim_res/<sim>/points
            candidates = [
                base / t,
                base / (t + ".pos"),
                self.proj.root / "sim_res" / sim / t,
                self.proj.root / "sim_res" / sim / (t + ".pos")
            ]
            hit = next((c for c in candidates if c.exists()), None)
            resolved.append(str(hit if hit else (base / (t + ".pos"))))

        cmd = [self.py(), self.runsie(), "post", sim, "-th", str(self.th_post.value())]
        if self.acc_post.isChecked():
            cmd += ["-a"]
        if self.etm_post.value() != 0:
            cmd += ["--etm", str(self.etm_post.value())]
        jobs_spec = self.jobs_post.text().strip()
        if self.sel_by_lambda_post.isChecked():
            lam_csv = self.lambda_post_edit.text().strip()
            if lam_csv:
                cmd += ["--lambdas", lam_csv]
        elif jobs_spec:
            cmd += ["-j", jobs_spec]
        cmd += ["-p"] + resolved

        self.runner.env.update(self._collect_env())
        self.runner.run(cmd)

    def do_points(self):
        sim = self.cur_sim()
        if not sim:
            QtWidgets.QMessageBox.warning(self, "No simulation", "Pick or create a simulation first.")
            return
        plane = self.plane.currentText()
        # compact any whitespace inside comma-separated lists (e.g. "-500, 500" -> "-500,500")
        def _csv(s: str) -> str:
            return re.sub(r"\s+", "", s.strip())
        args = [self.py(), self.runsie(), "make-points", sim, plane]
        # helper: add a closed interval only if both ends are provided
        def add_range(tag, mn, mx):
            nonlocal args
            tmin = mn.text().strip(); tmax = mx.text().strip()
            if tmin != "" and tmax != "":
                args += [f"--{tag}", tmin, tmax]
        # ranges
        add_range("x", self.xmin, self.xmax)
        add_range("y", self.ymin, self.ymax)
        add_range("z", self.zmin, self.zmax)
        # steps
        if self.step.text().strip(): args += ["--step", self.step.text().strip()]
        if self.stepx.text().strip(): args += ["--stepx", self.stepx.text().strip()]
        if self.stepy.text().strip(): args += ["--stepy", self.stepy.text().strip()]
        if self.stepz.text().strip(): args += ["--stepz", self.stepz.text().strip()]
        # fixed coordinates
        if plane == "xz" and self.y0.text().strip():
            args += [f"--y0={_csv(self.y0.text())}"]
        elif plane == "xy" and self.z0.text().strip():
            args += [f"--z0={_csv(self.z0.text())}"]
        elif plane == "yz" and self.x0.text().strip():
            args += [f"--x0={_csv(self.x0.text())}"]

        # Output name: add ".pos" if missing; if bare name, place under sim_res/<sim>/points/
        if self.outname.text().strip():
            name = self.outname.text().strip()
            # add extension if needed
            if not name.lower().endswith(".pos"):
                name += ".pos"
            # if no path given, save under sim_res/<sim>/points/
            p = Path(name)
            if not p.is_absolute() and (p.parent == Path(".")):
                base = self.proj.root / "sim_res" / sim / "points"
                base.mkdir(parents=True, exist_ok=True)
                name = str(base / p.name)
            args += ["-o", name]
        self.runner.run(args)

    def do_table(self):
        filename = self.tab_filename.text().strip()
        if not filename:
            QtWidgets.QMessageBox.warning(self, "No filename", "Please provide an output filename (e.g., erfc.bin).")
            return

        # Build the command using the CLI structure in run_sie.py
        cmd = [
            self.py(),
            self.runsie(),
            "make-table",
            filename,
            "--maxRe", str(self.tab_maxRe.value()),
            "--maxIm", str(self.tab_maxIm.value()),
            "--incRe", str(self.tab_incRe.value()),
            "--incIm", str(self.tab_incIm.value())
        ]

        if self.tab_text_mode.isChecked():
            cmd.append("--text")

        # Run the command in the background
        self.runner.env.update(self._collect_env())
        self.runner.run(cmd)

    def do_visualize(self):
        sim = self.cur_sim()
        if not sim:
            QtWidgets.QMessageBox.warning(self, "No simulation", "Pick or create a simulation first.")
            return
        # Decide save vs GUI-only preview
        persist = self.v_save.isChecked() and self.v_savepath.text().strip() != ""
        if persist:
            name = self.v_savepath.text().strip()
            # add .png if missing
            if not name.lower().endswith(".png"):
                name += ".png"
            # create media/ subfolder in sim_res/<sim>/
            media_dir = self.proj.root / "sim_res" / sim / "out" / "media"
            media_dir.mkdir(parents=True, exist_ok=True)
            out_png = media_dir / name
        else:
            out_png = self.proj.root / "__tmp_vis.png"

        # Build base command (note: visualization.py can accept multiple sims for --plot-csc)
        args = [self.py(), self.vispy()]

        # Branch on the selected modes
        if self.v_sel_xs.isChecked():
            # Cross sections (optionally with Mie overlay if parameters are filled)
            sims = [i.text() for i in self.v_simlist.selectedItems()] if hasattr(self, "v_simlist") else []
            if not sims:
                sims = [sim]
            args += sims + ["--root", "sim_res", "--plot-csc"]
            if self.v_subjob.value() > 0: args += ["--subjob", str(self.v_subjob.value())]
            # extract keys scs/acs/ext from labels
            sel_items = [i.text() for i in self.v_series.selectedItems()]
            series_keys = []
            for t in (sel_items or []):
                m = re.search(r"\((scs|acs|ext)\)", t, re.IGNORECASE)
                if m: series_keys.append(m.group(1).lower())
            if not series_keys:
                series_keys = ["scs", "acs", "ext"]
            for s in series_keys:
                args += ["--csc", s]
            # Mie overlay only if enabled AND user provides material + radius
            if self.v_mie_enable.isChecked() and self.v_mie_R.value() > 0 and self.v_mie_mat.text().strip():
                args += ["--mie-radius", f"{self.v_mie_R.value():.9g}",
                         "--mie-material", shlex.quote(self.v_mie_mat.text().strip()),
                         "--mie-background", shlex.quote(self.v_mie_bg.text().strip() or "eps=1,0"),
                         "--materials-root", "materials"]
                if self.spline_interp.isChecked():
                    args += ["--mie-spline"]
        elif self.v_sel_rt.isChecked():
            # R/T vs λ
            args += [sim, "--root", "sim_res", "--plot-rt"]
            if self.v_rt_points.text().strip():
                pts = self.v_rt_points.text().strip()
                if not pts.lower().endswith(".pos"):
                    pts += ".pos"
                args += ["--points", pts]
            args += ["--subjob", str(self.v_subjob.value())]
            prop_sel = self.v_rt_prop.currentText()
            if prop_sel == "auto":
                args += ["--rt-prop", "auto"]
            else:
                args += ["--rt-prop", "+1" if prop_sel.startswith("+") else "-1"]
        else:
            # Near-fields
            if self.v_points.text().strip():
                args += [sim, "--root", "sim_res"]
                pts = self.v_points.text().strip()
                if not pts.lower().endswith(".pos"):
                    pts += ".pos"
                args += ["--points", pts]
            if self.v_field.currentText(): args += ["--field", self.v_field.currentText()]
            if self.v_part.currentText():  args += ["--part", self.v_part.currentText()]
            if self.v_quan.currentText():  args += ["--quantity", self.v_quan.currentText()]
            if self.v_func.currentText():  args += ["--function", self.v_func.currentText()]
            if self.v_mode.currentText():  args += ["--mode", self.v_mode.currentText()]
            if self.v_subjob.value() > 0:  args += ["--subjob", str(self.v_subjob.value())]
            if self.v_log.isChecked(): args += ["--log"]
            if self.v_line.text().strip(): args += ["--line", self.v_line.text().strip()]
            # Only add hlines if checked
            if self.v_show_ints.isChecked():
                # Auto-detect interface z-lines from config.txt (ok if absent)
                try:
                    root_res   = self.proj.root / "sim_res"
                    interfaces = resolve_interfaces(sim, root_res, None)
                    if interfaces:
                        val = ",".join(f"{z:.9g}" for z in interfaces)
                        args += [f"--hlines={val}"]
                except Exception:
                    pass
            if self.v_lambda.value() > 0:
                args += ["--lambda", f"{self.v_lambda.value():.9g}"]

        # Always render to a PNG for preview; delete afterwards if not persisting
        args += ["--save", str(out_png)]

        self.status.showMessage("Rendering via visualization.py …")
        self.log.appendPlainText("$ " + " ".join(args))
        ret = subprocess.run(" ".join(args), cwd=str(self.proj.root), shell=True).returncode
        self.log.appendPlainText(f"[visualization.py exit {ret}]")

        # Display the saved image in the right panel
        self.v_plot.clear()
        try:
            if not out_png.exists():
                raise FileNotFoundError(str(out_png))
            img = mpimg.imread(str(out_png))
            self.v_plot.ax.imshow(img)
            self.v_plot.ax.axis('off')
        except Exception as e:
            self.v_plot.ax.text(0.5, 0.5, f"No image produced. {e}", ha='center', va='center')
        self.v_plot.draw()

        if not persist:
            try:
                out_png.unlink(missing_ok=True)
            except Exception:
                pass

    def do_mesh(self):
        sim = self.cur_sim()
        if not sim:
            QtWidgets.QMessageBox.warning(self, "No simulation", "Pick or create a simulation first.")
            return

        # Resolve mesh + interfaces the same way the CLI does
        try:
            root_res   = self.proj.root / "sim_res"    # same default as CLI
            mesh_path  = auto_find_mesh(sim, root_res, None)
            interfaces = resolve_interfaces(sim, root_res, None)
            nodes, faces, doms = read_mesh_mesh(mesh_path)
        except Exception as e:
            QtWidgets.QMessageBox.warning(self, "Load error", str(e))
            return

        mode     = self.m_mode.currentText()
        decim    = max(1, int(self.m_decimate.value()))
        alpha    = float(self.m_alpha.value())

        # Switch axes projection to match mode
        self.m_fig.clear()
        if mode.startswith("3D"):
            ax = self.m_fig.add_subplot(111, projection="3d")
            self.m_ax = ax
            title = f"Mesh 3D — {sim}"
            # Use the same plotter the CLI uses (draws into our ax)
            try:
                plot_mesh_3d(nodes, faces, doms,
                            interfaces=interfaces, title=title,
                            decimate=decim, alpha_mesh=alpha, ax=ax)
            except Exception as e:
                QtWidgets.QMessageBox.warning(self, "Plot error (3D)", str(e))
                return

        elif mode.startswith("2D"):
            ax = self.m_fig.add_subplot(111)  # 2D axes
            self.m_ax = ax
            plane = mode.split(":")[1].strip()  # "XY" / "XZ" / "YZ"
            title = f"Mesh 2D — {sim}"
            try:
                self._mesh_draw_2d(ax, nodes, faces, doms, plane=plane, decimate=decim,
                                   interfaces=interfaces, title=title)
            except Exception as e:
                QtWidgets.QMessageBox.warning(self, "Plot error (2D)", str(e))
                return

        else:  # 1D: z-interfaces
            ax = self.m_fig.add_subplot(111)
            self.m_ax = ax
            try:
                self._mesh_draw_1d(ax, nodes, interfaces, title=f"Z distribution — {sim}")
            except Exception as e:
                QtWidgets.QMessageBox.warning(self, "Plot error (1D)", str(e))
                return

        # Optional save
        if self.m_save.isChecked() and self.m_savepath.text().strip():
            name = self.m_savepath.text().strip()
            # add .png if missing
            if not name.lower().endswith(".png"):
                name += ".png"
            # create media/ subfolder in sim_res/<sim>/out/
            media_dir = self.proj.root / "sim_res" / sim / "out" / "media"
            media_dir.mkdir(parents=True, exist_ok=True)
            out_png = media_dir / name
            try:
                self.m_fig.savefig(str(out_png), dpi=300, bbox_inches="tight")
                self.log.appendPlainText(f"[mesh] saved: {out_png}")
            except Exception as e:
                self.log.appendPlainText(f"[mesh] save error: {e}")        
        
        self.m_canvas.draw_idle()

    # ------------------------- mesh helpers -------------------------
    def _mesh_draw_2d(self, ax, nodes, faces, doms, *, plane: str, decimate: int, interfaces, title: str):
        """
        Fast 2D rendering of triangle edges projected on a principal plane.
        Colors/legend follow the same (front->back) grouping as 3D.
        """
        plane = plane.upper()  # "XY"/"XZ"/"YZ"
        # Choose coordinate indices for projection
        if plane == "XY":
            i, j = 0, 1   # plot x vs y, ignore z
            xlabel, ylabel = "x", "y"
        elif plane == "XZ":
            i, j = 0, 2
            xlabel, ylabel = "x", "z"
        else:  # "YZ"
            i, j = 1, 2
            xlabel, ylabel = "y", "z"

        ax.clear()
        ax.set_title(title)
        ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)

        # Group faces by transition and build line segments per group
        groups = build_transition_groups(faces, doms)
        keys   = sorted(groups.keys())
        cols   = color_wheel(len(keys))
        handles, labels = [], []

        for ci, key in enumerate(keys):
            idxs = decimate_indices(groups[key], decimate)
            if idxs.size == 0:
                continue
            tri = faces[idxs]  # (K,3) node indices

            # Build segments for edges of each triangle in 2D projection
            # each tri -> 3 segments: (a,b), (b,c), (c,a)
            A = nodes[tri[:,0]]; B = nodes[tri[:,1]]; C = nodes[tri[:,2]]
            segs = []
            segs.extend(np.stack([A[:,[i,j]], B[:,[i,j]]], axis=1))
            segs.extend(np.stack([B[:,[i,j]], C[:,[i,j]]], axis=1))
            segs.extend(np.stack([C[:,[i,j]], A[:,[i,j]]], axis=1))
            segs = np.asarray(segs)

            lc = LineCollection(segs, colors=[cols[ci]], linewidths=0.25, alpha=0.9)
            ax.add_collection(lc)

            handles.append(Patch(facecolor=cols[ci], edgecolor="k", alpha=0.7))
            labels.append(f"{key[0]}→{key[1]}")

        ax.autoscale()
        ax.set_aspect("equal", adjustable="box")
        ax.grid(True, alpha=0.25)

        # Draw interface reference lines (if provided)
        if interfaces:
            for zi in interfaces:
                ax.axhline(y=zi, color="0.35", lw=0.6, ls="--", alpha=0.8)

        if handles:
            leg = ax.legend(handles=handles, labels=labels,
                            title="Domain boundaries (front→back):\nbackground domain = 0",
                            loc="upper left", frameon=True)
            try:
                leg._legend_box.align = "left"
            except Exception:
                pass

    def _mesh_draw_1d(self, ax, nodes, interfaces, title: str):
        """
        1D view: histogram of node z distribution with horizontal lines at interface z values.
        (Simple, informative; matches layered-media focus.)
        """
        ax.clear()
        ax.set_title(title)
        z = nodes[:,2]
        ax.hist(z, bins=50, density=False, alpha=0.6)
        ax.set_xlabel("z")
        ax.set_ylabel("node count")
        ax.grid(True, alpha=0.3)
        if interfaces:
            for zi in interfaces:
                ax.axvline(x=zi, color="0.35", lw=0.8, ls="--", alpha=0.9)

    def do_distributed(self):
        sim = self.cur_sim()
        if not sim:
            QtWidgets.QMessageBox.warning(self, "No simulation", "Pick or create a simulation first.")
            return
        # map GUI mode to numeric --mode for the Python distributed runner
        mode_text = self.dist_mode.currentText()
        if mode_text == "isolated scatterer(s) in homogeneous background":
            mode_idx = 0
        elif mode_text == "periodic structure in homogeneous background":
            mode_idx = 1
        else:
            mode_idx = 2

        runner_path = self.dist_runner.text().strip()
        if not runner_path:
            QtWidgets.QMessageBox.warning(
                self,
                "Missing distributed runner",
                "Provide the path to dist_run_helios.py in Settings."
            )
            return

        root = str(self.proj.root)
        # Call the Python-based runner
        cmd = [
            self.py(),
            runner_path,
            "--mode", str(mode_idx),
            "-s", sim,
            "--root", root,
            "-R", self.run_sie_edit.text().strip(),
        ]

        nodes_csv = self.dist_nodes.text().strip()
        if not nodes_csv:
            QtWidgets.QMessageBox.warning(self, "No nodes", "Enter node names or enable 'Use count'.")
            return
        cmd += ['-N', nodes_csv]

        sub = self.dist_subcmd.currentText()
        if sub == 'Run solver':
            cmd += ['solve', '-l', str(self.dist_level.value())]
        elif sub == 'Run post-processing':
            if not self.dist_points.text().strip():
                QtWidgets.QMessageBox.warning(
                    self,
                    "No points",
                    "Provide one or more point-set names for distributed post-processing."
                )
                return
            cmd += ['post']
            pts_tokens = self.dist_points.text().strip().split()
            resolved = []
            base = Path(root).expanduser() / 'sim_res' / sim / 'points'
            for t in pts_tokens:
                tpath = Path(t)
                if tpath.is_absolute() and tpath.exists():
                    resolved.append(str(tpath))
                    continue
                candidates = [
                    base / t,
                    base / (t + ".pos"),
                    base.parent / t,
                    base.parent / (t + ".pos")
                ]
                hit = next((c for c in candidates if c.exists()), None)
                resolved.append(str(hit if hit else (base / (t + ".pos"))))
            cmd += ['-p'] + resolved

        if mode_idx == 1:  # periodic homogeneous
            cmd += ["--etm", str(self.dist_etm.value())]
        if self.dist_sel_by_lambda.isChecked() and self.dist_lambda_edit.text().strip():
            cmd += ["--lambdas", self.dist_lambda_edit.text().strip()]
        elif self.dist_jobs.text().strip():
            cmd += ["-j", self.dist_jobs.text().strip()]
        if self.dist_acc.isChecked(): cmd += ['-a']
        cmd += ['-th', str(self.dist_threads.value())]

        self.runner.env.update(self._collect_env())
        self.runner.run(cmd)

    def _tab_settings(self):
        w = QtWidgets.QWidget(); f = QtWidgets.QFormLayout(w)
        self.py_exec = QtWidgets.QLineEdit(sys.executable)
        self.run_sie_edit = QtWidgets.QLineEdit(str(self.proj.run_sie))
        self.vis_edit = QtWidgets.QLineEdit(str(self.proj.vis))
        # Optional shell launchers used by Distributed tab
        self.dist_runner = QtWidgets.QLineEdit(str(self.proj.root / "dist_run_helios.py"))
        self.env_table = QtWidgets.QTableWidget(0, 2); self.env_table.setHorizontalHeaderLabels(["KEY","VALUE"])
        self.env_table.horizontalHeader().setStretchLastSection(True)
        self.btn_add_env = QtWidgets.QPushButton("Add environment variable")
        self.btn_add_env.clicked.connect(self.add_env_row)
        f.addRow("Python:", self.py_exec)
        f.addRow("run_sie.py:", self.run_sie_edit)
        f.addRow("visualization.py:", self.vis_edit)
        f.addRow("Distributed runner (dist_run_helios.py):", self.dist_runner)
        f.addRow(self.btn_add_env)
        f.addRow(self.env_table)
        return w

    def _bind_lambda_toggle(self, check: QtWidgets.QCheckBox, jobs_row: QtWidgets.QWidget, lam_row: QtWidgets.QWidget):
        """Show 'jobs_row' when unchecked, 'lam_row' when checked; initialize state."""
        def _apply():
            use_l = check.isChecked()
            jobs_row.setVisible(not use_l)
            lam_row.setVisible(use_l)
        check.toggled.connect(_apply)
        _apply()

# ----------------------------- entry -----------------------------
def main():
    ap = argparse.ArgumentParser(description="HELIOS GUI")
    ap.add_argument("--root", default=".", help="Project root (folder containing run_sie.py, pytools/, sim_data/, sim_res/)")
    args = ap.parse_args()

    app = QtWidgets.QApplication(sys.argv)
    proj = ProjectPaths.discover(Path(args.root))
    logo = Path(proj.root).resolve() / "logo.png"
    icon = QtGui.QIcon(str(logo))
    if not icon.isNull():
        app.setWindowIcon(icon)          # application icon (dock/alt-tab)
    else:
        print("icon not found:", logo)

    # then create the window
    gui = HeliosGUI(proj)

    # mirror app icon to the window (helps on some WMs)
    if not app.windowIcon().isNull():
        gui.setWindowIcon(app.windowIcon())

    gui.show()
    sys.exit(app.exec())

if __name__ == "__main__":
    main()
