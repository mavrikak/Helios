#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_sie.py — Orchestrate HELIOS workflows.

This script provides a subcommand-based CLI that prepares simulations,
runs the C++ solvers, and post-processes fields for visualization. It is
designed to work with the HELIOS project layout and with C++11 builds of 
SIENano / SIENanoPP.

───────────────────────────────────────────────────────────────────────────────
Project layout (relative to repository root)
───────────────────────────────────────────────────────────────────────────────
sim_data/<sim>/       # user inputs (config.txt, mesh.mphtxt)
sim_res/<sim>/        # generated outputs for this simulation
  ├─ jobs/            # job.N.xml files (one per wavelength)
  ├─ mesh.mesh        # converted mesh (from COMSOL *.mphtxt or Gmsh *.msh)
  ├─ out/             # solver outputs
  │   ├─ csc/         # normalized location of *.csc files (gathered here)
  │   └─ fields/      # field files produced by SIENanoPP (.ein,.esc,.hin,.hsc)
  ├─ points/          # sampling point sets (*.pos) produced by make-points
  └─ logs/            # one log per solver run

Optional materials/tables:
- materials/          # ε(λ) tables used by jobwriter
- erfc.bin            # periodic erfc lookup table (auto-resolved if needed)

───────────────────────────────────────────────────────────────────────────────
Subcommands (high level)
───────────────────────────────────────────────────────────────────────────────
prepare     Generate jobs/ from config.txt and convert mesh if needed.
solve       Run the SIENano app over selected jobs.
post        Run the SIENanoPP app on produced *.sol using user-supplied points.
make-points Create regular sampling grids; avoids layers' interfaces.
all         prepare -> solve (-> post if --points/-p is given).

───────────────────────────────────────────────────────────────────────────────
Typical workflow
───────────────────────────────────────────────────────────────────────────────
1) Put inputs in sim_data/<sim>/:
     - config.txt          (see project examples)
     - Mesh file:
        * COMSOL text mesh (*.mphtxt), or
        * Gmsh MSH v4 ASCII (*.msh).
NOTE: COMSOL and Gmsh inputs support 3D tetrahedral meshes only.

2) Prepare jobs and mesh:
     run_sie.py prepare <sim> [--mode 0|1|2] [--meshfile 0|1]

3) Solve:
     run_sie.py solve   <sim> [-a] [-l LEVEL] [-th THREADS] [--jobs 1,3,5-8]
                        [--table PATH]

4) Post-process fields:
     run_sie.py post    <sim> -p points/<file.pos> [-a] [-th THREADS] 
                        [--table PATH]

───────────────────────────────────────────────────────────────────────────────
Periodic simulation awareness
───────────────────────────────────────────────────────────────────────────────
- Jobs are treated as periodic if their job XML contains a <periodic> block.
- For periodic jobs, this script resolves an erfc lookup table and passes both:
    * environment: ERFC_TABLE=/abs/path/to/erfc.bin
    * CLI flag:    -t /abs/path/to/erfc.bin
  You may override with --table PATH. If unset, common locations are tried
  under sim_res/<sim>/ and the project root.

───────────────────────────────────────────────────────────────────────────────
make-points behavior (interface-safe sampling)
───────────────────────────────────────────────────────────────────────────────
- Generates plane grids (xz, xy, yz) with independent steps per axis.
- Loads interface z-levels from sim_res/<sim>/config.txt (if present) and
  nudges any point that would land exactly on an interface by 1e-2 (display
  units unchanged). This avoids some degenerate evaluations in SIENanoPP.

───────────────────────────────────────────────────────────────────────────────
I/O conventions
───────────────────────────────────────────────────────────────────────────────
- Points files (*.pos): plain text, one "x y z" per line (no header).
- Solutions (*.sol): produced by SIENano and kept directly under
  sim_res/<sim>/out/ so SIENanoPP can consume them without shims.
- Cross sections (*.csc): produced by SIENano; this script moves them into
  sim_res/<sim>/out/csc/ so visualization can find them.
- Field files: SIENanoPP writes .ein/.esc/.hin/.hsc into out/; this script
  organizes them under out/fields/<points_name>/ per input points set.
───────────────────────────────────────────────────────────────────────────────
CLI synopsis (selected options)
───────────────────────────────────────────────────────────────────────────────
Common:
  --root ROOT            Project root (auto-detected in most setups)
  --launcher CMD         Prefix external runner, e.g. "srun -N 1 -n 16"
  -a, --accurate         Enable solver accuracy options (passes -a)
  -th, --threads N       Threads for SIENano/SIENanoPP
  --table PATH           erfc table path for periodic jobs
  --jobs SPEC            Job selection: "1,3,5-8" (1-based)

prepare:
  --mode {0,1,2}         0: isolated homogeneous, 1: periodic homogeneous,
                         2: isolated layered
  --meshfile {0,1}       1 converts the unique *.mphtxt or *.msh into *.mesh

make-points:
  --plane {xz,xy,yz}             Target plane
  --step H                       Uniform step (overridden by per-axis)
  --stepx/--stepy/--stepz        Per-axis steps
  --x0 --x1 --y0 --y1 --z0 --z1  Bounds in native units
  -p, --points NAME              Output basename under sim_res/<sim>/points/

post:
  -p, --points PATH|NAME Points file or name under sim_res/<sim>/points/

───────────────────────────────────────────────────────────────────────────────
Environment variables
───────────────────────────────────────────────────────────────────────────────
NO_BAR=1                 Disable progress bars/spinners (CI-friendly)
ERFC_TABLE=/path/file    Default periodic table (overridden by --table)

───────────────────────────────────────────────────────────────────────────────
Exit codes & logging
───────────────────────────────────────────────────────────────────────────────
- Non-zero exit codes indicate user-facing errors (bad paths/flags/inputs).
- All external tool invocations are logged under sim_res/<sim>/logs/.
- On non-interactive TTYs, spinners are suppressed to keep logs clean.

"""

from __future__ import annotations
import argparse
import os
import re
import shutil
import sys, time
from tqdm import tqdm
from pathlib import Path
from typing import List, Optional, Tuple
import subprocess
import shlex

# ------------------------------- paths ---------------------------------------

ROOT    = Path(__file__).resolve().parent
APPS    = ROOT / "apps"
PYTOOLS = ROOT / "pytools"

SIENANO   = APPS / "SIENano"
SIENANOPP = APPS / "SIENanoPP"
SIETABLE  = APPS / "SIETable"

# ----------------------------- utilities -------------------------------------

def die(msg: str, code: int = 2) -> None:
    """
    Terminate the program with a short, prefixed error message and given exit code.
    Use for user-facing, expected failure modes (bad paths/flags/missing files).
    """
    print(f"Error: {msg}", file=sys.stderr)
    sys.exit(code)

def ensure_exists(p: Path, kind: str = "path") -> None:
    """Guard helper: exit with a clear message if a required file/folder is missing."""
    if not p.exists():
        die(f"Cannot find {kind}: {p}")

def list_jobs(jobs_dir: Path) -> List[Path]:
    """Return job.*.xml files sorted by numeric index (job.<N>.xml -> N asc)."""
    return sorted(
        (p for p in jobs_dir.glob("job.*.xml") if p.is_file()),
        key=lambda p: int(p.stem.split(".")[1])
    )

def rel_or_abs(p: Path, base: Path) -> Path:
    """Interpret 'p' relative to 'base' if not absolute; return a Path."""
    return p if p.is_absolute() else (base / p)

def find_unique_meshfile(sim_data_dir: Path) -> Optional[Path]:
    """
    Prefer a single mesh source under sim_data/<sim>: one of
      - *.msh (Gmsh, v4 ASCII recommended) or
      - *.mphtxt (COMSOL)
    If multiple are present, error out for explicit user choice.
    """
    cands = sorted(list(sim_data_dir.glob("*.msh")) + list(sim_data_dir.glob("*.mphtxt")))
    if not cands:
        return None
    if len(cands) > 1:
        die(f"Found multiple mesh sources in {sim_data_dir}. Use --meshfile to disambiguate.")
    return cands[0]

def tee_run(cmd: List[str], cwd: Path, log_file: Path, show_prefix: str = "") -> int:
    """
    Run a command, tee stdout/stderr to a log file, and (optionally) show a
    minimal progress indication on TTY. On non-TTY (e.g., ssh/CI), we suppress
    spinners to avoid noisy output.

    Notes:
      - Appends to the log file (does not truncate). Each call writes a '$ cmd' header.
      - Respect SIE_NO_SPINNER=1 to disable the lightweight TTY prefix line.
    """
    log_file.parent.mkdir(parents=True, exist_ok=True)
    is_tty = sys.stdout.isatty()
    no_spinner = os.environ.get("SIE_NO_SPINNER", "0") == "1"

    # Append instead of write: keeps any header callers added before tee_run()
    with log_file.open("a") as lf:
        lf.write(f"$ {' '.join(cmd)}\n\n")
        lf.flush()

        proc = subprocess.Popen(
            cmd, cwd=str(cwd),
            stdout=lf, stderr=subprocess.STDOUT,
            universal_newlines=True
        )

        # Non-interactive context (CI/ssh) or user disabled spinner -> no TTY noise.
        if not is_tty or no_spinner:
            # Just wait quietly; avoid carriage returns/garbage in non-TTY
            ret = proc.wait()
            if is_tty and show_prefix:
                # one clean line after completion
                print(f"{show_prefix} done (rc={ret})")
            return ret

        # Light TTY indicator (no animated spinner flood)
        if show_prefix:
            # Minimal on-tty progress (single-line prefix) instead of an animated spinner.
            print(f"{show_prefix} ...", end="", flush=True)

        ret = proc.wait()

        if show_prefix:
            print(f"\r{show_prefix} done (rc={ret}){' ' * 10}")
        return ret

_bars = {}
def progress_bar(i: int, n: int, prefix: str = "", width: int = 40) -> None:
    """
    tqdm-based progress bar that:
      - shows in interactive terminals
      - renders on stderr (so stdout prints remain intact)
      - persists at the bottom like tqdm, even while logs scroll
    Keeps the same call signature as before: progress_bar(i, n, prefix).
    """
    global _bars
    is_tty = sys.stdout.isatty()
    no_bar = bool(os.environ.get("NO_BAR", False))

    # Fallback to a simple completion line when not interactive or disabled
    if no_bar or not is_tty or n <= 1:
        if i == n:
            print(f"{prefix}{i}/{n} done.")
        return

    key = prefix.strip() or "default"
    bar = _bars.get(key)
    if bar is None:
        # Leave=True keeps the final bar on screen; dynamic_ncols adapts to terminal width.
        bar = tqdm(
            total=n,
            desc=key,
            dynamic_ncols=True,
            leave=True,
            file=sys.stderr,
            bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt}",
        )
        _bars[key] = bar

    # Update only the delta from current count
    delta = max(0, i - bar.n)
    if delta:
        bar.update(delta)

    if i >= n:
        bar.close()
        _bars.pop(key, None)
        # ensure a clean line after the closed bar
        print("", file=sys.stderr)

def parse_launcher(launcher: str | None) -> list[str]:
    """
    Parse an optional launcher prefix (e.g., 'srun -N 1 -n 16').
    Uses shlex.split to respect quotes. Validates the binary exists on PATH.
    """
    if not launcher:
        return []
    parts = shlex.split(launcher)
    if not parts:
        return []
    exe = parts[0]
    if shutil.which(exe) is None:
        raise FileNotFoundError(
            f"Launcher executable not found: '{exe}'. "
            f"Install it, load its module, provide an absolute path, "
            f"or drop --launcher to run locally."
        )
    return parts

# ---- points (.pos) helpers ---------------------------------------------------

def _frange(a: float, b: float, h: float, eps: float = 1e-9) -> list[float]:
    """
    Inclusive float range [a, b] with fixed step h.
    Mitigates FP drift by rounding step count and tolerating a small epsilon near 'b'.
    """
    if h <= 0:
        die(f"step must be > 0, got {h}")
    
    # Round step count so tiny FP discrepancies don't drop the last point.
    n = int(round((b - a) / h)) + 1
    vals = []
    for i in range(max(0, n)):
        v = a + i * h
        if v > b + eps:
            break
        vals.append(v)
    # If 'b' is still meaningfully beyond the last value, force-append it.
    if vals and (b - vals[-1]) > (h * 0.5) and (b - vals[-1]) > eps:
        vals.append(b)
    return vals

def write_pos_grid(sim: str, plane: str,
                   x: tuple[float, float] | None,
                   y: tuple[float, float] | None,
                   z: tuple[float, float] | None,
                   x0: float | None, y0: float | None, z0: float | None,
                   out_name: str,
                   *,                           # keyword-only extras
                   stepx: float | None = None,
                   stepy: float | None = None,
                   stepz: float | None = None,
                   interfaces: set[float] = frozenset(),
                   append: bool = False) -> None:
    """
    Generate a regular grid of points on a named plane and save to
    sim_res/<sim>/points/<out_name>.pos as plain text:
        x y z
        x y z
        ...
    """
    res_dir = ROOT / "sim_res" / sim
    ensure_exists(res_dir, "sim_res/<sim> folder")
    pts_dir = res_dir / "points"
    pts_dir.mkdir(parents=True, exist_ok=True)
    out_path = pts_dir / out_name
    mesh_path = res_dir / "mesh.mesh"
    z_limits = get_mesh_z_limits(mesh_path)

    pts: list[tuple[float, float, float]] = []

    if plane == "xz":
        if x is None or z is None or y0 is None:
            die("xz requires --x xmin xmax --z zmin zmax --y0 value")
        xs = _frange(x[0], x[1], stepx)
        zs = _frange(z[0], z[1], stepz)
        # Load nodes on this y-plane (axis 1)
        bad_2d = load_planar_nodes(mesh_path, 1, float(y0))
        bad_xs = {p[0] for p in bad_2d}
        bad_zs = {p[1] for p in bad_2d}
        x_mid = (x[0] + x[1]) * 0.5
        z_mid = (z[0] + z[1]) * 0.5
        for xi in xs:
            if round(xi, 4) in bad_xs:
                if xi >= x_mid: xi -= 0.5 * stepx
                else: xi += 0.5 * stepx
            for zi in zs:
                zp = nudge_z_if_interface(zi, interfaces, z_limits, 0.5 * stepz)
                pts.append((xi, float(y0), zp))

    elif plane == "xy":
        if x is None or y is None or z0 is None:
            die("xy requires --x xmin xmax --y ymin ymax --z0 value")
        xs = _frange(x[0], x[1], stepx)
        ys = _frange(y[0], y[1], stepy)
        # Load nodes on this z-plane (axis 2)
        bad_2d = load_planar_nodes(mesh_path, 2, float(z0))
        bad_xs = {p[0] for p in bad_2d}
        bad_ys = {p[1] for p in bad_2d}
        x_mid = (x[0] + x[1]) * 0.5
        y_mid = (y[0] + y[1]) * 0.5
        z_fixed = nudge_z_if_interface(float(z0), interfaces, z_limits)
        if abs(z_fixed) < 1e-2: z_fixed = z_fixed + 1e-2
        for xi in xs:
            if round(xi, 4) in bad_xs:
                if xi >= x_mid: xi -= 0.5 * stepx
                else: xi += 0.5 * stepx
            for yi in ys:
                if round(yi, 4) in bad_ys:
                    if yi >= y_mid: yi -= 0.5 * stepy
                    else: yi += 0.5 * stepy
                pts.append((xi, yi, z_fixed))

    elif plane == "yz":
        if y is None or z is None or x0 is None:
            die("yz requires --y ymin ymax --z zmin zmax --x0 value")
        ys = _frange(y[0], y[1], stepy)
        zs = _frange(z[0], z[1], stepz)
        # Load nodes on this x-plane (axis 0)
        bad_2d = load_planar_nodes(mesh_path, 0, float(x0))
        bad_ys = {p[0] for p in bad_2d}
        bad_zs = {p[1] for p in bad_2d}
        y_mid = (y[0] + y[1]) * 0.5
        z_mid = (z[0] + z[1]) * 0.5
        for yi in ys:
            if round(yi, 4) in bad_ys:
                if yi >= y_mid: yi -= 0.5 * stepy
                else: yi += 0.5 * stepy
            for zi in zs:
                zp = nudge_z_if_interface(zi, interfaces, z_limits, 0.5 * stepz)
                pts.append((float(x0), yi, zp))
    else:
        die(f"Unknown plane '{plane}' (choose from xz, xy, yz)")

    with out_path.open("a" if append else "w") as f:
        for xp, yp, zp in pts:
            f.write(f"{xp:.6f} {yp:.6f} {zp:.6f}\n")

    print(f"[points] Wrote {len(pts)} points -> {out_path}")

def parse_job_selection(spec: str, max_job: int) -> list[int]:
    """
    Parse job selection strings (e.g., '1,3,5-8') into sorted, 1-based indices.
    Clips ranges to [1, max_job] and rejects empty/invalid specs.
    """
    sel: set[int] = set()
    for tok in spec.split(","):
        tok = tok.strip()
        if not tok:
            continue
        # Range token (A-B): accept reversed bounds and clip to valid job interval.
        if "-" in tok:
            a, b = tok.split("-", 1)
            a, b = int(a), int(b)
            lo, hi = min(a, b), max(a, b)
            for k in range(max(1, lo), min(max_job, hi) + 1):
                sel.add(k)
        else:
            k = int(tok)
            if 1 <= k <= max_job:
                sel.add(k)
    out = sorted(sel)
    if not out:
        die(f"No valid jobs in selection '{spec}' (1..{max_job}).")
    return out

# ---- solution file helpers (per job, with subjobs) --------------------------

SOL_RE = re.compile(r"^sim(\d+)_inc_(\d+)\.sol$")

def list_job_solutions(sol_dir: Path, job_number: int) -> list[Path]:
    """
    Return all solutions for a given job_number as a list of Paths, sorted by subjob index.
    Only files matching sim<i>_inc_<j>.sol are considered, and only those with i==job_number.
    """
    found: list[tuple[int, Path]] = []
    if not sol_dir.exists():
        return []
    for p in sol_dir.glob("*.sol"):
        m = SOL_RE.match(p.name)
        if not m:
            continue
        i = int(m.group(1))
        sj = int(m.group(2))
        if i == job_number:
            found.append((sj, p))
    return [p for _, p in sorted(found)]

# -------------------------- periodic helpers ---------------------------------

PERIODIC_TAG_RE = re.compile(r"<\s*periodic\b", re.IGNORECASE)

def job_is_periodic(job_xml: Path) -> bool:
    """
    Heuristic: treat a job as periodic if '<periodic' appears in its XML (case-insensitive).
    """
    try:
        txt = job_xml.read_text(errors="ignore")
    except Exception:
        return False
    return PERIODIC_TAG_RE.search(txt) is not None

def extract_periodic_params_from_config(cfg: Path) -> Optional[List[float]]:
    """
    Scan config.txt for a line with five numerics: px py cx cy zcut (commas/spaces allowed).
    Returns [px, py, cx, cy, zcut] as floats if found, else None.
    """
    params = None
    # use UTF-8 with BOM support and ignore stray bytes
    with cfg.open("r", encoding="utf-8-sig", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = re.split(r"[,\s]+", s)
            if len(parts) == 5:
                try:
                    vals = list(map(float, parts))
                except ValueError:
                    continue
                params = vals
    return params

# ---- table resolution (works with or without --table) -----------------------

def resolve_table_path_optional(sim: str, res_dir: Path, table_opt: Optional[str]) -> Optional[str]:
    """
    Resolve a periodic table path to an absolute path, or None.
    Priority:
      1) --table value, interpreted as:
         - absolute path, or
         - relative to sim_res/<sim>, or
         - relative to project ROOT, or
         - relative to current working directory
      2) ERFC_TABLE environment variable (absolute or any of the above)
      3) sim_res/<sim>/erfc.bin
      4) ROOT/third_party/erfc.bin
      5) ROOT/erfc.bin
    Returns an absolute string path if found; otherwise None.
    """
    candidates: list[Path] = []

    # Treat user paths as either absolute or relative to sim_res/<sim>, project ROOT, or CWD.
    def add_candidates_from(s: Optional[str]):
        if not s:
            return
        p = Path(s)
        if p.is_absolute():
            candidates.append(p)
        else:
            candidates.extend([res_dir / p, ROOT / p, Path.cwd() / p])

    # Common placements used in this project (local per-sim or shared under project tree):
    # 1) CLI
    add_candidates_from(table_opt)
    # 2) ENV
    add_candidates_from(os.environ.get("ERFC_TABLE"))
    # 3) Common local placements
    candidates.append(res_dir / "erfc.bin")
    candidates.append(ROOT / "third_party" / "erfc.bin")
    candidates.append(ROOT / "erfc.bin")

    for c in candidates:
        try:
            if c.exists() and c.is_file():
                return str(c.resolve())
        except Exception:
            continue
    return None

def ensure_table_for_periodic(sim: str, res_dir: Path, table_opt: Optional[str]) -> str:
    """
    Ensure we have a valid table path for periodic jobs. If found, export ERFC_TABLE
    and return absolute path. If not found, error.
    """
    tab_abs = resolve_table_path_optional(sim, res_dir, table_opt)
    if not tab_abs:
        die(
            "Periodic job detected but no periodic lookup table was found.\n"
            "Provide one via:\n"
            "  --table /abs/path/to/erfc.bin   (or path relative to project or sim_res/<sim>)\n"
            "or set the environment variable:\n"
            "  export ERFC_TABLE=/abs/path/to/erfc.bin"
        )
    # Export for child processes (SIENano/SIENanoPP) even if user didn't set it.
    os.environ["ERFC_TABLE"] = tab_abs  # make it visible to children
    return tab_abs

# --- Interfaces parsing from sim_res/<sim>/config.txt ------------------------

def load_interface_z(res_dir: Path) -> set[float]:
    """
    Parse the 'Layered media interfaces' section in sim_res/<sim>/config.txt.
    Collect numeric z-values on subsequent non-empty, non-comment lines 
    (stops on blank/comment).
    """
    zs: set[float] = set()
    cfg = res_dir / "config.txt"
    if not cfg.exists():
        return zs
    lines = cfg.read_text(errors="ignore").splitlines()
    for i, line in enumerate(lines):
        # Start collecting immediately after this marker line; 
        # accept comma- or space-separated numbers.
        if "Layered media interfaces" in line:
            j = i + 2
            while j < len(lines):
                s = lines[j].strip()
                if not s or s.startswith("#"):
                    break
                for tok in re.split(r"[,\s]+", s):
                    if not tok:
                        continue
                    try:
                        zs.add(float(tok))
                    except ValueError:
                        pass
                j += 1
            break
    return zs

def nudge_z_if_interface(z: float, interfaces: set[float], 
                         z_limits: tuple[float, float] | None = None, 
                         offset: float = 1e-2,
                         tol: float = 1e-6) -> float:
    """
    If z coincides with an interface (within 'tol'), shift it to avoid singular sampling.
    Logic:
      - If <= z_mid (mesh mid): shift up     (+1e-2)
      - If >= z_mid (mesh mid): shift down   (-1e-2)
      - Otherwise:              shift down   (-1e-2)
    """
    for zi in interfaces:
        if z_limits and abs(z - zi) < tol:
            z_min, z_max = z_limits
            z_mid = (z_min + z_max) * 0.5
            if zi <= z_mid:
                return z + offset
            # Default to negative shift
            return z - offset
    return z

def load_planar_nodes(mesh_path: Path, plane_idx: int, plane_val: float, tol: float = 1e-4) -> set[tuple[float, float]]:
    """
    Reads mesh.mesh to find nodes lying on the specific plane (axis[plane_idx] ~= plane_val).
    Returns a set of (u, v) coordinates (rounded to 4 decimals) for the other two axes.
    Assumes simple mesh format: lines with exactly 3 floats are nodes.
    """
    bad_pts = set()
    if not mesh_path.exists():
        return bad_pts
    try:
        with mesh_path.open("r") as f:
            for line in f:
                parts = line.split()
                if len(parts) == 3:  # Heuristic: 3 columns = x y z
                    try:
                        c = [float(p) for p in parts]
                        if abs(c[plane_idx] - plane_val) < tol:
                            # Extract the other two coordinates in index order
                            others = [c[i] for i in range(3) if i != plane_idx]
                            bad_pts.add((round(others[0], 4), round(others[1], 4)))
                    except ValueError:
                        pass
    except Exception:
        pass
    return bad_pts
    
def get_mesh_z_limits(mesh_path: Path) -> tuple[float, float] | None:
    """
    Scan mesh.mesh for global min/max Z coordinates (3rd column).
    Returns (z_min, z_max) or None if failed.
    """
    if not mesh_path.exists():
        return None
    zmin, zmax = float('inf'), -float('inf')
    try:
        with mesh_path.open("r") as f:
            for line in f:
                parts = line.split()
                if len(parts) == 3:
                    try:
                        z = float(parts[2])
                        if z < zmin: zmin = z
                        if z > zmax: zmax = z
                    except ValueError:
                        pass
    except Exception:
        return None
    if zmin == float('inf'):
        return None
    return zmin, zmax

# ---- wavelength helpers (lambdalist.txt and job.N.xml parsing) ---------------
LAMBDA_TAG_RE = re.compile(r"<\s*(wavelength|lambda|lambda0)[^>]*>([^<]+)<", re.IGNORECASE)
NUM_UNIT_RE = re.compile(r"([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s*(nm|um|µm|micron|m|ev|hz|khz|mhz|ghz|thz)?", re.IGNORECASE)

def parse_lambda_list(expr: str) -> List[float]:
    """
    Parse a wavelength list expression into sorted unique nm values.
    Accepted forms (comma/semicolon separated):
      690
      690nm
      0.69um
      690-696
      690:2:696     (start:step:end)
    Tokens can mix units; all are converted to nm.
    Returns a sorted list of unique wavelengths in nm.
    """
    if not expr:
        return []
    toks = [t.strip() for t in expr.replace(";", ",").split(",") if t.strip()]
    vals: List[float] = []

    def parse_one(s: str) -> Optional[float]:
        m = NUM_UNIT_RE.fullmatch(s)
        if not m:
            return None
        v = float(m.group(1))
        u = (m.group(2) or "nm")
        return _to_nm(v, u)

    for t in toks:
        # step range a:s:b ?
        if ":" in t:
            parts = [p.strip() for p in t.split(":")]
            if len(parts) == 3:
                a = _to_nm(*NUM_UNIT_RE.fullmatch(parts[0]).groups()) if NUM_UNIT_RE.fullmatch(parts[0]) else None
                s = _to_nm(*NUM_UNIT_RE.fullmatch(parts[1]).groups()) if NUM_UNIT_RE.fullmatch(parts[1]) else None
                b = _to_nm(*NUM_UNIT_RE.fullmatch(parts[2]).groups()) if NUM_UNIT_RE.fullmatch(parts[2]) else None
                # step given in wavelength units (nm after conversion)
                if a is not None and s is not None and b is not None and s != 0:
                    cur = a
                    # handle direction via sign of step
                    step = s if b >= a else -abs(s)
                    # include end if it lands within epsilon
                    eps = abs(step) * 1e-9 + 1e-9
                    if step > 0:
                        while cur <= b + eps:
                            vals.append(cur)
                            cur += step
                    else:
                        while cur >= b - eps:
                            vals.append(cur)
                            cur += step
                    continue
        # simple range a-b ?
        if "-" in t and not t.strip().startswith("-"):
            a_str, b_str = [x.strip() for x in t.split("-", 1)]
            a = parse_one(a_str)
            b = parse_one(b_str)
            if a is not None and b is not None:
                step = 1.0  # default 1 nm
                if b >= a:
                    cur = a
                    while cur <= b + 1e-9:
                        vals.append(cur)
                        cur += step
                else:
                    cur = a
                    while cur >= b - 1e-9:
                        vals.append(cur)
                        cur -= step
                continue
        # single value
        v = parse_one(t)
        if v is not None:
            vals.append(v)

    # normalize: sort + unique (within 1e-6 nm tolerance)
    vals.sort()
    out: List[float] = []
    for v in vals:
        if not out or abs(out[-1] - v) > 1e-6:
            out.append(v)
    return out

def _to_nm(val: float, unit: str | None) -> float | None:
    if not unit or unit.lower() == "nm":
        return val
    u = unit.lower()
    if u in ("um", "µm", "micron"): return val * 1e3
    if u == "m": return val * 1e9
    if u == "ev": return 1239.841984 / val if val != 0 else None
    fact = {"hz":1, "khz":1e3, "mhz":1e6, "ghz":1e9, "thz":1e12}.get(u)
    if fact:
        f = val * fact
        c = 299_792_458.0
        return (c / f) * 1e9 if f != 0 else None
    return None

def _parse_lambda_from_xml(xml: Path) -> float | None:
    try:
        txt = xml.read_text(errors="ignore")
    except Exception:
        return None
    m = LAMBDA_TAG_RE.search(txt)
    def first_nm(s: str) -> float | None:
        m2 = NUM_UNIT_RE.search(s)
        if not m2: return None
        val, unit = float(m2.group(1)), (m2.group(2) or "nm")
        return _to_nm(val, unit)
    if m:
        nm = first_nm(m.group(2))
        if nm is not None: return nm
    m3 = NUM_UNIT_RE.search(txt)
    if m3:
        return _to_nm(float(m3.group(1)), m3.group(2) or "nm")
    return None

def ensure_lambdalist(res_dir: Path) -> dict[int, float]:
    """Return {job_index -> λ_nm}. Writes lambdalist.txt if missing."""
    lf = res_dir / "lambdalist.txt"
    idx2wl: dict[int, float] = {}
    if lf.exists():
        for i, line in enumerate(lf.read_text().splitlines(), 1):
            s = line.strip()
            if not s or s[0] in "#%/": continue
            parts = s.split()
            if len(parts) == 1:
                idx2wl[i] = float(parts[0])
            elif len(parts) >= 2:
                # 'job.N  λ' or similar
                try:
                    n = int(parts[0].split(".")[-1])
                    idx2wl[n] = float(parts[1])
                except Exception:
                    pass
        if idx2wl:
            return idx2wl

    # Build from jobs/*.xml
    jobs_dir = res_dir / "jobs"
    jobs = list_jobs(jobs_dir)
    if not jobs:
        die(f"No jobs in {jobs_dir} to infer wavelengths from.")
    for p in jobs:
        n = int(p.stem.split(".")[1])  # job.N.xml -> N
        nm = _parse_lambda_from_xml(p)
        if nm is not None:
            idx2wl[n] = nm
    if not idx2wl:
        die("Could not infer wavelengths from job XMLs.")
    with lf.open("w") as f:
        for n in sorted(idx2wl):
            f.write(f"job.{n} {idx2wl[n]:.9g}\n")
    return idx2wl

def lambdas_to_job_indices(res_dir: Path, lam_expr: str) -> List[int]:
    """Map lambda expression (supports numbers, a-b, a:s:b, mixed units) to nearest job indices."""
    idx2wl = ensure_lambdalist(res_dir)
    if not idx2wl:
        die("No wavelengths found to map --lambdas.")
    idxs  = sorted(idx2wl.keys())
    wls   = [idx2wl[i] for i in idxs]

    asks = parse_lambda_list(lam_expr)
    if not asks:
        die("Empty --lambdas list.")
    chosen: set[int] = set()
    for lam in asks:
        j = min(range(len(idxs)), key=lambda k: abs(wls[k] - lam))
        chosen.add(idxs[j])
    return sorted(chosen)

# ------------------------ prepare logic (with prompt) ------------------------

def _confirm_delete(res_dir: Path) -> bool:
    """Interactive safeguard before deleting an existing sim_res/<sim> folder."""
    print(f"Folder already exists: {res_dir}")
    ans = input("Delete it and proceed? [y/N]: ").strip().lower()
    return ans in ("y", "yes")

def _maybe_reset_res(sim: str, overwrite: bool, keep: bool) -> Path:
    """Create/refresh sim_res/<sim> depending on --overwrite/--keep and TTY interactivity."""
    res_dir = ROOT / "sim_res" / sim
    if not res_dir.exists():
        return res_dir
    if overwrite:
        shutil.rmtree(res_dir)
        return res_dir
    if keep:
        return res_dir
    # if non-interactive (stdin not a TTY), default to keep to avoid EOF
    try:
        import sys
        # In non-interactive mode (e.g., batch), default to KEEP to avoid EOF on input().
        if not sys.stdin.isatty():
            return res_dir
    except Exception:
        return res_dir
    if _confirm_delete(res_dir):
        shutil.rmtree(res_dir)
        return res_dir
    else:
        die("Aborted by user (kept existing sim_res).")

def prepare_sim(sim: str, mode: int, meshfile: Optional[str], overwrite: bool, keep: bool, 
                materials: str = "eps", spline: bool = False) -> Tuple[Path, Path]:
    """
    Ensure sim_res/<sim> structure exists and is populated:
    - Copy config.txt for provenance.
    - If jobs/ is empty, generate jobs via pytools/jobwriter.py (mode, meshfile, sim).
    - If meshfile==1 and mesh.mesh is missing, convert the lone *.mphtxt via pytools/meshconvert.py.
    When mode==1 (periodic) and px py cx cy zcut exist in config, pass them to meshconvert.
    """
    sim_data = ROOT / "sim_data" / sim
    ensure_exists(sim_data, "sim_data folder")
    cfg = sim_data / "config.txt"
    ensure_exists(cfg, "config.txt")

    res_dir = _maybe_reset_res(sim, overwrite, keep)
    jobs_dir = res_dir / "jobs"
    out_dir  = res_dir / "out"
    logs_dir = res_dir / "logs"
    for d in (res_dir, jobs_dir, out_dir, logs_dir):
        d.mkdir(parents=True, exist_ok=True)
    
    # Keep a snapshot of the config next to results to make runs self-contained/reproducible.
    cfg_dest = res_dir / "config.txt"
    try:
        # copy if missing or source is newer
        shutil.copy2(str(cfg), str(cfg_dest))
        print(f"[prepare] Copied config.txt -> {cfg_dest}")
    except Exception as e:
        print(f"[prepare] Warning: could not copy config.txt: {e}")

    # 1) JOBS via jobwriter.py (if empty)
    jobs = list_jobs(jobs_dir)
    if not jobs:
        print(f"[prepare] Generating jobs from {cfg} ...")
        jobwriter = PYTOOLS / "jobwriter.py"
        ensure_exists(jobwriter, "jobwriter.py")
        cmd = [sys.executable, str(jobwriter), str(mode), "0", sim]
        if spline: cmd.append("--spline")
        env = os.environ.copy()
        if materials == "nk":
            env["MATERIALS_COLUMNS"] = "nk"
        else:
            # ensure explicit eps mode (unset to let jobwriter default to eps)
            env.pop("MATERIALS_COLUMNS", None)
        ret = subprocess.run(cmd, cwd=str(PYTOOLS), env=env).returncode
        if ret != 0:
            die("jobwriter failed.")
        jobs = list_jobs(jobs_dir)
        if not jobs:
            die("jobwriter finished but no jobs were produced.")

    # 2) MESH: reuse existing *.mesh if present, otherwise call meshconvert.py
    mesh_out = res_dir / "mesh.mesh"
    if not mesh_out.exists():
        # If a .mesh already exists in sim_data/<sim>, just move it
        existing_meshes = sorted(sim_data.glob("*.mesh"))
        if existing_meshes:
            if len(existing_meshes) > 1:
                die(f"Multiple .mesh files found in {sim_data}; not sure which to use.")
            src = existing_meshes[0]
            print(f"[prepare] Found existing mesh {src.name}; moving to {mesh_out.name} ...")
            shutil.copy(str(src), str(mesh_out))
            print(f"[prepare] Done. jobs: {len(jobs)}  mesh: {mesh_out.exists()}")
            return res_dir, jobs_dir

        # otherwise fall back to converting a .msh / .mphtxt via meshconvert.py
        meshconvert = PYTOOLS / "meshconvert.py"
        ensure_exists(meshconvert, "meshconvert.py")
        # Resolve mesh source
        src = None
        if meshfile:
            # accept absolute or relative to sim_data/<sim> or project root
            cand = Path(meshfile)
            if not cand.is_absolute():
                for base in (sim_data, ROOT):
                    p = base / meshfile
                    if p.exists():
                        cand = p; break
            if not cand.exists():
                die(f"--meshfile not found: {meshfile}")
            src = cand
        else:
            src = find_unique_meshfile(sim_data)
            if src is None:
                die(f"No mesh file found in {sim_data}. Provide --meshfile or place a *.msh / *.mphtxt there.")

        # Copy the source mesh next to results
        src_dest = res_dir / src.name
        if src_dest.resolve() != src.resolve():
            try: shutil.copy2(str(src), str(src_dest))
            except Exception: shutil.copyfile(str(src), str(src_dest))

        # Apply periodic cropping (mode==1) for both COMSOL and Gmsh sources
        if mode == 1:
            pp = extract_periodic_params_from_config(cfg)
            if pp and len(pp) == 5:
                print(f"[prepare] Converting mesh with periodic crop {pp}: {src_dest.name} -> mesh.mesh ...")
                cmd = [sys.executable, str(meshconvert), "-o", str(mesh_out), str(src_dest)] + list(map(str, pp))
            else:
                print("[prepare] Periodic mode but params not found in config; converting as isolated.")
                cmd = [sys.executable, str(meshconvert), "-o", str(mesh_out), str(src_dest)]
        else:
            print(f"[prepare] Converting mesh: {src_dest.name} -> mesh.mesh ...")
            cmd = [sys.executable, str(meshconvert), "-o", str(mesh_out), str(src_dest)]

        ret = subprocess.run(cmd, cwd=str(PYTOOLS)).returncode
        if ret != 0:
            die("meshconvert failed.")

    print(f"[prepare] Done. jobs: {len(jobs)}  mesh: {mesh_out.exists()}")
    return res_dir, jobs_dir

# --------------------------------- solve -------------------------------------

def run_solve(sim: str, mode: int, accurate: bool, 
              level: int, threads: int, launcher: Optional[str],
              etm: int, table_opt: Optional[str],
              jobs_spec: Optional[str] = None) -> None:
    """
    Run SIENano over selected jobs (or all) under sim_res/<sim>/jobs/.
    - Adds -a/-l/-th as requested.
    - Detects periodic jobs via XML and appends -etm and -t <table>.
    - Moves produced sim<i>_inc_<j>.sol into out/sol/ and keeps logs per job.
    """
    res_dir = ROOT / "sim_res" / sim
    jobs_dir = res_dir / "jobs"
    out_dir = res_dir / "out"
    logs_dir = res_dir / "logs"
    out_dir.mkdir(exist_ok=True, parents=True)
    logs_dir.mkdir(exist_ok=True, parents=True)

    jobs_all = list_jobs(jobs_dir)
    if not jobs_all:
        die(f"No jobs in {jobs_dir}.")

    total_job_count = len(jobs_all)
    lamcsv = getattr(sys.modules[__name__], "_args_lambdas", None)

    def _safe_parse(spec: Optional[str]) -> List[int]:
        if not spec:
            return []
        try:
            return parse_job_selection(spec, max_job=total_job_count)  # 1-based indices
        except SystemExit:
            return []
        except Exception:
            return []

    # λ -> indices (sorted unique)
    mapped_idxs: Optional[List[int]] = lambdas_to_job_indices(res_dir, lamcsv) if lamcsv else None
    # slice indices from launcher (-j)
    slice_idxs: Optional[List[int]] = _safe_parse(jobs_spec) if jobs_spec else None

    def _is_contiguous_block(idx: List[int]) -> bool:
        if not idx:
            return False
        return sorted(idx) == list(range(min(idx), min(idx) + len(idx)))

    # Decide final indices
    if mapped_idxs is not None and slice_idxs is not None:
        # Prefer round-robin distribution across slices when the slice is a contiguous block.
        if _is_contiguous_block(slice_idxs):
            block_size = len(slice_idxs)
            if block_size <= 0:
                print("[solve] Empty slice; exiting.", flush=True)
                return
            group_id = (min(slice_idxs) - 1) // block_size
            total_groups = max(1, (total_job_count + block_size - 1) // block_size)
            # Round-robin over the λ-filtered list so slices split evenly
            final_idxs = mapped_idxs[group_id::total_groups]
            if not final_idxs:
                print("[solve] No jobs for this slice after lambda mapping; exiting.", flush=True)
                return
        else:
            # Fallback: intersect by numeric ids
            final_idxs = sorted(set(mapped_idxs).intersection(slice_idxs))
            if not final_idxs:
                print("[solve] No jobs for this slice after lambda mapping; exiting.", flush=True)
                return
    elif mapped_idxs is not None:
        final_idxs = mapped_idxs
    elif slice_idxs is not None:
        if not slice_idxs:
            print("[solve] No valid jobs in slice; exiting.", flush=True)
            return
        final_idxs = slice_idxs
    else:
        final_idxs = list(range(1, total_job_count + 1))

    if lamcsv: print(f"[solve] lambda selection mapped to jobs: {final_idxs}", flush=True)

    # materialize job paths from indices
    jobs: List[Path] = [(res_dir / "jobs" / f"job.{i}.xml") for i in final_idxs]

    base_flags: List[str] = []
    if accurate:
        base_flags += ["-a"]
    base_flags += ["-l", str(level)]
    base_flags += ["-th", str(threads)]

    launch_prefix = parse_launcher(launcher)

    n = len(jobs)
    print(f"[solve] Running {n} job(s) for '{sim}'\n")
    print(f"[solve-progress] total {n}", flush=True)
    for idx, job in enumerate(jobs, 1):
        job_number = int(job.stem.split(".")[1])  # from job.<N>.xml
        progress_bar(idx - 1, n, prefix="[solve] ")

        # per-job periodic flags
        per_flags: List[str] = []
        if job_is_periodic(job):
            per_flags += ["-etm", str(etm)]
            tab_abs = ensure_table_for_periodic(sim, res_dir, table_opt)
            per_flags += ["-t", tab_abs]

        cmd = launch_prefix + [str(SIENANO), "-j", str(job.relative_to(res_dir))] + base_flags + per_flags
        log_name = f"S{job_number}.txt"
        ret = tee_run(cmd, cwd=res_dir, log_file=logs_dir / log_name,
                      show_prefix=f"[solve job {idx}/{n}]")
        if ret != 0:
            die(f"SIENano failed for {job.name}. See logs/{log_name}")

        # Collect only this job's cross-section files into out/csc/
        csc_dir = out_dir / "csc"
        csc_dir.mkdir(exist_ok=True, parents=True)
        # Cross sections are named sim<i>_inc_<j>.csc where i == job_number
        for csc in out_dir.glob(f"sim{job_number}_inc_*.csc"):
            dest = csc_dir / csc.name
            if dest.resolve() != csc.resolve():
                try:
                    shutil.move(str(csc), str(dest))
                except Exception:
                    # If another node moved it first, ignore
                    pass
        progress_bar(idx, n, prefix="[solve] ")
        print(f"[solve-progress] done {idx}", flush=True)
    
    print("\n[solve] Done.")

# --------------------------------- post --------------------------------------

def run_post(sim: str, pos_files: List[Path],
             accurate: bool, threads: int,
             launcher: Optional[str], etm: int, table_opt: Optional[str],
             jobs_spec: Optional[str] = None) -> None:
    """
    Run SIENanoPP for each (selected job) x (provided .pos set).
    - Matches sim<i>_inc_<j>.sol by job number i (all subjobs j).
    - Writes logs per (job, points) and moves new field files into out/fields/<points>/.
    - Periodic-aware: adds -etm/-t when the originating job was periodic.
    - Safe for concurrent nodes: only the delta of new files is moved after each run.
    """
    res_dir = ROOT / "sim_res" / sim
    ensure_exists(res_dir, "sim_res/<sim> folder")
    jobs_dir = res_dir / "jobs"
    out_dir  = res_dir / "out"
    logs_dir = res_dir / "logs"
    ensure_exists(jobs_dir, "jobs folder")
    ensure_exists(out_dir,  "out folder")
    logs_dir.mkdir(exist_ok=True, parents=True)

    jobs_all = list_jobs(jobs_dir)
    if not jobs_all:
        die(f"No jobs in {jobs_dir}.")

    total_job_count = len(jobs_all)
    lamcsv = getattr(sys.modules[__name__], "_args_lambdas", None)

    def _safe_parse(spec: Optional[str]) -> List[int]:
        if not spec:
            return []
        try:
            return parse_job_selection(spec, max_job=total_job_count)
        except SystemExit:
            return []
        except Exception:
            return []

    mapped_idxs: Optional[List[int]] = lambdas_to_job_indices(res_dir, lamcsv) if lamcsv else None
    slice_idxs:  Optional[List[int]] = _safe_parse(jobs_spec) if jobs_spec else None

    def _is_contiguous_block(idx: List[int]) -> bool:
        if not idx:
            return False
        return sorted(idx) == list(range(min(idx), min(idx) + len(idx)))

    if mapped_idxs is not None and slice_idxs is not None:
        if _is_contiguous_block(slice_idxs):
            block_size = len(slice_idxs)
            if block_size <= 0:
                print("[post] Empty slice; exiting.", flush=True)
                return
            group_id = (min(slice_idxs) - 1) // block_size
            total_groups = max(1, (total_job_count + block_size - 1) // block_size)
            final_idxs = mapped_idxs[group_id::total_groups]
            if not final_idxs:
                print("[post] No jobs for this slice after lambda mapping; exiting.", flush=True)
                return
        else:
            final_idxs = sorted(set(mapped_idxs).intersection(slice_idxs))
            if not final_idxs:
                print("[post] No jobs for this slice after lambda mapping; exiting.", flush=True)
                return
    elif mapped_idxs is not None:
        final_idxs = mapped_idxs
    elif slice_idxs is not None:
        if not slice_idxs:
            print("[post] No valid jobs in slice; exiting.", flush=True)
            return
        final_idxs = slice_idxs
    else:
        final_idxs = list(range(1, total_job_count + 1))

    jobs: List[Path] = [(res_dir / "jobs" / f"job.{i}.xml") for i in final_idxs]
    if not jobs:
        print("[post] No jobs selected; exiting.", flush=True)
        return

    # Resolve .pos paths
    pos_files_resolved = [rel_or_abs(p, res_dir) for p in pos_files]
    for pth in pos_files_resolved:
        ensure_exists(pth, ".pos file")

    # Solutions now live directly under out/, so we search there.
    sol_dir = out_dir
    total = len(jobs) * len(pos_files_resolved) # One SIENanoPP call per (job x points)
    if total == 0:
        die(f"No matching solutions sim<i>_inc_<j>.sol under {sol_dir} for the selected job(s).")

    # Base flags and launcher
    base_flags: List[str] = []
    if accurate:
        base_flags += ["-a"]
    base_flags += ["-th", str(threads)]
    launch_prefix = parse_launcher(launcher)

    k = 0
    print(f"[post] {len(jobs)} job(s) x {len(pos_files_resolved)} point-set(s) = {total} runs\n")
    print(f"[post-progress] total {total}", flush=True)

    for job in jobs:
        job_number = int(job.stem.split(".")[1])  # from "job.<N>.xml"

        # Periodic flags per job
        per_flags: List[str] = []
        if job_is_periodic(job):
            per_flags += ["-etm", str(etm)]
            tab_abs = ensure_table_for_periodic(sim, res_dir, table_opt)
            per_flags += ["-t", tab_abs]

        # Find any one solution file for this job
        sols_for_job = list_job_solutions(sol_dir, job_number)
        if not sols_for_job:
            print(f"[post] Warning: no .sol for job {job_number} in {sol_dir}")
            continue
        first_sol = sols_for_job[0]

        for pos in pos_files_resolved:
            points_name = pos.stem

            # Per-(job, points) log file
            log_name = f"P{job_number}_{points_name}.txt"
            log_path = logs_dir / log_name
            with log_path.open("a") as lf:
                lf.write(f"# Post-processing logs for {job.name} @ points '{points_name}'\n\n")

            # Output fields directory for this points set
            fields_dir = out_dir / "fields" / points_name
            fields_dir.mkdir(parents=True, exist_ok=True)

            k += 1
            progress_bar(k - 1, total, prefix="[post] ")

            sol_shim = out_dir / first_sol.name

            before = {p.name for p in out_dir.glob("*")}

            cmd = launch_prefix + [
                str(SIENANOPP),
                "-j", str(job.relative_to(res_dir)),
                "-s", str(sol_shim.relative_to(res_dir)),
                "-p", str(pos.relative_to(res_dir))
            ] + base_flags + per_flags

            ret = tee_run(cmd, cwd=res_dir, log_file=log_path,
                            show_prefix=f"[post job {job_number} x {points_name}]")
            if ret != 0:
                die(f"SIENanoPP failed for {job.name} with points '{points_name}'. See logs/{log_name}")

            after = {p.name for p in out_dir.glob("*")}
            new_names = after - before
            exts = (".ein", ".esc", ".hin", ".hsc")
            for name in sorted(new_names):
                if not name.endswith(exts):
                    continue
                f = out_dir / name
                dest = fields_dir / f.name
                if dest.resolve() == f.resolve():
                    continue
                try:
                    shutil.move(str(f), str(dest))
                except Exception:
                    pass

            print(f"[post-progress] done {k}", flush=True)
            progress_bar(k, total, prefix="[post] ")

    print("\n[post] Done (one SIENanoPP call per job and point-set).")

def run_make_table(filename: str, as_text: bool, maxRe: float, maxIm: float, incRe: float, incIm: float) -> None:
    """Run SIETable to generate a periodic lookup table."""
    ensure_exists(SIETABLE, "apps/SIETable")
    
    flag = "-w" if as_text else "-wb"
    cmd = [str(SIETABLE), flag, filename]
    
    # maintable.cpp expects these 4 values from stdin
    input_str = f"{maxRe} {maxIm} {incRe} {incIm}\n"
    
    print(f"[make-table] Generating table: {filename}")
    print(f"[make-table] Inputs: maxRe={maxRe}, maxIm={maxIm}, incRe={incRe}, incIm={incIm}")
    
    try:
        subprocess.run(
            cmd,
            input=input_str,
            text=True,
            check=True
        )
        
        # Move the file to third_party/
        source_path = Path(filename)
        if source_path.exists():
            dest_dir = ROOT / "third_party"
            dest_dir.mkdir(parents=True, exist_ok=True) # Ensure folder exists
            dest_path = dest_dir / source_path.name
            
            # Move the file (overwrites if it already exists there)
            shutil.move(str(source_path), str(dest_path))
            
            print(f"[make-table] Successfully created and moved to: {dest_path.relative_to(ROOT)}")
        else:
            print(f"[make-table] Warning: Expected output file {filename} not found.")
            
    except subprocess.CalledProcessError as e:
        die(f"SIETable failed with exit code {e.returncode}")

# ----------------------------------- CLI -------------------------------------

def main():
    p = argparse.ArgumentParser(
        description="Run jobs from a sim_data subfolder (periodic-ready).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    # Subcommand-based CLI to mirror the typical workflow (prepare -> solve -> post).
    sub = p.add_subparsers(dest="cmd", required=True)

    # Common flags shared by prepare/solve/post/all.
    def add_common_args(sp):
        sp.add_argument("sim", help="sim_data/<sim> subfolder name")
        sp.add_argument("--mode", type=int, default=2, choices=(0,1,2),
                        help="0: iso-hom | 1: per-hom | 2: iso-layered")
        sp.add_argument("--meshfile", type=str, default=None,
                        help="Path to input mesh (.msh v4 ASCII or .mphtxt). If omitted, auto-detect under sim_data/<sim>.")
        sp.add_argument("-a", "--accurate", action="store_true",
                        help="pass -a to SIENano / SIENanoPP")
        sp.add_argument("-l", "--level", type=int, default=2,
                        help="SIENano -l LEVEL (default 2)")
        sp.add_argument("-th", "--threads", type=int, default=0,
                        help="SIENano/SIENanoPP -th THREADS (default 0)")
        sp.add_argument("--launcher", type=str, default=None,
                        help='command prefix for nodes/queues, e.g. "srun -N 1 -n 16"')
        sp.add_argument("--overwrite", action="store_true",
                        help="if sim_res/<sim> exists, DELETE it without prompting")
        sp.add_argument("--keep", action="store_true",
                        help="if sim_res/<sim> exists, KEEP it without prompting")
        # periodic extras
        sp.add_argument("--etm", type=int, default=1,
                        help="periodic evaluation terms (forwarded as -etm) when job is periodic")
        sp.add_argument("--table", type=str, default=None,
                        help="path to periodic erfc lookup table; if omitted we auto-detect")

    # <subcommand>: brief summary of what the stage does and how it interacts with others.
    p_prep = sub.add_parser("prepare", help="Generate jobs and mesh (asks before deleting existing sim_res)")
    add_common_args(p_prep)
    p_prep.add_argument("--materials", choices=("eps", "nk"), default="eps",
                        help="materials table interpretation for jobwriter: eps=(Re,Im) or nk=(n,k)")
    p_prep.add_argument("--spline", action="store_true",
                        help="use cubic spline interpolation for material tables (jobwriter --spline)")                        

    # <subcommand>: brief summary of what the stage does and how it interacts with others.
    p_solve = sub.add_parser("solve", help="Run SIENano on all jobs (auto-prepare; asks before deleting sim_res)")
    add_common_args(p_solve)
    p_solve.add_argument("-j", "--jobs", type=str, default=None,
                         help='job selection like "1,3,5-8" (defaults to all)')
    p_solve.add_argument("--lambdas", type=str, default=None,
                         help='comma/semicolon-separated wavelengths (nm) to select closest matching jobs')                         

    # <subcommand>: brief summary of what the stage does and how it interacts with others.
    p_post = sub.add_parser("post", help="Run SIENanoPP on all jobs with given .pos files")
    add_common_args(p_post)
    p_post.add_argument("-p", "--pos", nargs="+", required=True,
                        help="one or more .pos files (relative to sim_res/<sim> or absolute)")
    p_post.add_argument("-j", "--jobs", type=str, default=None,
                        help='job selection like "1,3,5-8" (defaults to all)')
    p_post.add_argument("--lambdas", type=str, default=None,
                        help='comma/semicolon-separated wavelengths (nm) to select closest matching jobs')

    # <subcommand>: brief summary of what the stage does and how it interacts with others.
    p_all = sub.add_parser("all", help="prepare + solve (+ post if -p is given)")
    add_common_args(p_all)
    p_all.add_argument("--materials", choices=("eps", "nk"), default="eps",
                       help="materials table interpretation for jobwriter: eps=(Re,Im) or nk=(n,k)")
    p_all.add_argument("--spline", action="store_true",
                       help="use cubic spline interpolation for material tables (jobwriter --spline)")    
    p_all.add_argument("-p", "--pos", nargs="*", default=[],
                       help="optional .pos files to also run post-processing")
    p_all.add_argument("-j", "--jobs", type=str, default=None,
                       help='job selection like "1,3,5-8" applied to solve (and post if -p is given)')
    p_all.add_argument("--lambdas", type=str, default=None,
                       help='comma/semicolon-separated wavelengths (nm) to select closest matching jobs')
    
    # <subcommand>: brief summary of what the stage does and how it interacts with others.
    p_pts = sub.add_parser("make-points", help="Generate regular grid .pos files")
    p_pts.add_argument("sim", help="sim_data/<sim> subfolder name")
    p_pts.add_argument("plane", choices=("xz", "xy", "yz"), help="plane to sample")
    p_pts.add_argument("--x", nargs=2, type=float, metavar=("XMIN", "XMAX"), help="x-range")
    p_pts.add_argument("--y", nargs=2, type=float, metavar=("YMIN", "YMAX"), help="y-range")
    p_pts.add_argument("--z", nargs=2, type=float, metavar=("ZMIN", "ZMAX"), help="z-range")

    # step options
    p_pts.add_argument("--step",  type=float, help="uniform grid spacing for all axes")
    p_pts.add_argument("--stepx", type=float, help="grid spacing along x (overrides --step for x)")
    p_pts.add_argument("--stepy", type=float, help="grid spacing along y (overrides --step for y)")
    p_pts.add_argument("--stepz", type=float, help="grid spacing along z (overrides --step for z)")

    # Accept multiple values via comma/space separated list; parsed later
    p_pts.add_argument("--x0", type=str, help="fixed x value(s) for yz plane (comma/space separated)")
    p_pts.add_argument("--y0", type=str, help="fixed y value(s) for xz plane (comma/space separated)")
    p_pts.add_argument("--z0", type=str, help="fixed z value(s) for xy plane (comma/space separated)")
    p_pts.add_argument("-o", "--output", default=None, help="output filename (.pos)")

    # <subcommand>: generate lookup tables
    p_tab = sub.add_parser("make-table", help="Generate an erfc lookup table (binary or text)")
    p_tab.add_argument("filename", help="Output filename (e.g., erfc.bin)")
    p_tab.add_argument("--text", action="store_true", help="Write as a text table instead of binary (-w instead of -wb)")
    p_tab.add_argument("--maxRe", type=float, default=10.0, help="Maximum real part")
    p_tab.add_argument("--maxIm", type=float, default=10.0, help="Maximum imaginary part")
    p_tab.add_argument("--incRe", type=float, default=0.1, help="Increment for real part")
    p_tab.add_argument("--incIm", type=float, default=0.1, help="Increment for imaginary part")

    args = p.parse_args()

    # make --lambdas visible to run_solve/run_post without changing many signatures
    setattr(sys.modules[__name__], "_args_lambdas", getattr(args, "lambdas", None))

    # Sanity checks: binaries/scripts must exist relative to project ROOT.
    ensure_exists(SIENANO, "apps/SIENano")
    if not os.access(SIENANO, os.X_OK):
        die(f"apps/SIENano is not executable: {SIENANO}")    
    ensure_exists(SIENANOPP, "apps/SIENanoPP")
    if not os.access(SIENANOPP, os.X_OK):
        die(f"apps/SIENanoPP is not executable: {SIENANOPP}")
    ensure_exists(PYTOOLS / "jobwriter.py", "pytools/jobwriter.py")
    ensure_exists(PYTOOLS / "meshconvert.py", "pytools/meshconvert.py")

    if args.cmd == "prepare":
        prepare_sim(args.sim, args.mode, args.meshfile, args.overwrite, 
                    args.keep, materials=args.materials, spline=args.spline)

    elif args.cmd == "solve":
        run_solve(args.sim, args.mode, args.accurate, 
                  args.level, args.threads, args.launcher,
                  args.etm, args.table,
                  jobs_spec=args.jobs)

    elif args.cmd == "make-points":
        out_name = args.output or f"plane_{args.plane}.pos"
        xr = tuple(args.x) if args.x else None
        yr = tuple(args.y) if args.y else None
        zr = tuple(args.z) if args.z else None

        # Parse possibly multiple fixed-plane values (comma/space separated)
        def _parse_multi(txt: str | None) -> list[float]:
            if not txt:
                return []
            vals: list[float] = []
            for tok in re.split(r"[,\s]+", txt.strip()):
                if not tok:
                    continue
                try:
                    vals.append(float(tok))
                except ValueError:
                    die(f"Invalid numeric in fixed-plane list: '{tok}'")
            return vals

        # Derive per-axis step: a specific step overrides the uniform --step for that axis.
        def _axis_step(uniform, specific, axis):
            v = specific if specific is not None else uniform
            if v is None:
                die(f"--step{axis} (or --step) is required for {args.plane} plane")
            if v <= 0:
                die(f"--step{axis} must be > 0")
            return float(v)

        if args.plane == "xz":
            sx = _axis_step(args.step, args.stepx, "x")
            sz = _axis_step(args.step, args.stepz, "z")
            sy = None  # not used
            fixed_vals = _parse_multi(args.y0)
            fixed_tag  = "y0"
        elif args.plane == "xy":
            sx = _axis_step(args.step, args.stepx, "x")
            sy = _axis_step(args.step, args.stepy, "y")
            sz = None  # not used
            fixed_vals = _parse_multi(args.z0)
            fixed_tag  = "z0"
        else:  # yz
            sy = _axis_step(args.step, args.stepy, "y")
            sz = _axis_step(args.step, args.stepz, "z")
            sx = None  # not used
            fixed_vals = _parse_multi(args.x0)
            fixed_tag  = "x0"

        if not fixed_vals:
            die(f"{args.plane} requires --{fixed_tag} value(s).")

        # Load interface z-levels so generated grids avoid sampling exactly on layer boundaries.
        res_dir = ROOT / "sim_res" / args.sim
        interfaces = load_interface_z(res_dir)

        # Aggregate: write all requested planes into the SAME output file
        for i, v in enumerate(fixed_vals):
            x0v = v if fixed_tag == "x0" else None
            y0v = v if fixed_tag == "y0" else None
            z0v = v if fixed_tag == "z0" else None
            write_pos_grid(
                args.sim, args.plane,
                xr, yr, zr,
                x0=x0v, y0=y0v, z0=z0v,
                out_name=out_name,
                stepx=sx, stepy=sy, stepz=sz,
                interfaces=interfaces,
                append=(i > 0)
            )

    elif args.cmd == "post":
        pos_files = [Path(s) for s in args.pos]
        run_post(args.sim, pos_files, args.accurate, args.threads, args.launcher,
                 args.etm, args.table, jobs_spec=args.jobs)

    elif args.cmd == "all":
        prepare_sim(args.sim, args.mode, args.meshfile, args.overwrite, 
                    args.keep, materials=args.materials, spline=args.spline)
        run_solve(args.sim, args.mode, args.accurate, 
                  args.level, args.threads, args.launcher,
                  args.etm, args.table,
                  jobs_spec=args.jobs)
        if args.pos:
            pos_files = [Path(s) for s in args.pos]
            run_post(args.sim, pos_files, args.accurate, args.threads, args.launcher,
                     args.etm, args.table, jobs_spec=args.jobs)
    
    elif args.cmd == "make-table":
        run_make_table(args.filename, args.text, args.maxRe, args.maxIm, args.incRe, args.incIm)

if __name__ == "__main__":
    main()
