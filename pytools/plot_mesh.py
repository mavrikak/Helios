#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_mesh.py — 3D visualization of mesh.mesh + layered-media interfaces (if existent).

Usage examples
--------------
# From project layout (auto: sim_res/<sim>/mesh.mesh, interfaces from sim_res/<sim>/config.txt)
python3 pytools/plot_mesh.py --sim new_sim --root sim_res --show

# Direct file and explicit interfaces (nm units)
python3 pytools/plot_mesh.py --mesh path/to/mesh.mesh --interfaces 0,50,100 --show

# Save a PNG without opening a window
MPLBACKEND=Agg python3 pytools/plot_mesh.py --sim new_sim --save mesh_view.png
"""
from __future__ import annotations
import argparse
from pathlib import Path
from typing import List, Optional, Tuple, Dict

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.patches import Patch

# --------------------------- I/O helpers ---------------------------

def read_mesh_mesh(path: Path) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Read HELIOS mesh format:
      Nnodes
      x y z           (per node)
      Nfaces
      i j k front back   (1-based indices)

    Returns
    -------
    nodes : (N,3) float64
    faces : (M,3) int32   (0-based indices of nodes)
    doms  : (M,2) int32   (front, back)
    """
    with path.open("r", encoding="utf-8") as f:
        # nodes
        n_nodes = int(f.readline().strip())
        nodes = []
        for _ in range(n_nodes):
            x, y, z = map(float, f.readline().split())
            nodes.append((x, y, z))
        nodes = np.asarray(nodes, dtype=float)

        # faces
        n_faces = int(f.readline().strip())
        faces = []
        doms = []
        for _ in range(n_faces):
            i, j, k, fr, bk = f.readline().split()
            faces.append((int(i) - 1, int(j) - 1, int(k) - 1))  # 0-based
            doms.append((int(fr), int(bk)))
        faces = np.asarray(faces, dtype=np.int32)
        doms = np.asarray(doms, dtype=np.int32)

    return nodes, faces, doms


def parse_interfaces_from_config(config_path: Path) -> Optional[List[float]]:
    """
    Read 'Layered media interfaces' (bottom to top) from the config file.

    Rules:
      - Find a line containing 'Layered media interfaces'.
      - Then skip blank/comment/separator lines until a value line.
      - Parse numbers from one or more subsequent lines (commas or spaces).
      - Allow inline comments (anything after '#').
      - After values have started, stop when encountering a blank line
        or a comment line (i.e., next section header).
    """
    if not config_path.exists():
        return None

    vals: List[float] = []
    seen_header = False       # saw the section title
    started_values = False    # consumed at least one numeric token

    try:
        with config_path.open("r", encoding="utf-8") as f:
            for raw in f:
                line = raw.strip()

                # Look for the target header
                if not seen_header and "Layered media interfaces" in line:
                    seen_header = True
                    continue

                if not seen_header:
                    continue

                # Handle inline comments within the section
                if "#" in line:
                    line = line.split("#", 1)[0].strip()

                # Before values start: keep skipping empties/separators
                if not started_values and not line:
                    continue

                # If we've started collecting, an empty line ends the section
                if started_values and not line:
                    break

                # Normalize delimiters and parse tokens
                tokens = line.replace(",", " ").split()
                parsed_any = False
                for tok in tokens:
                    try:
                        vals.append(float(tok))
                        parsed_any = True
                    except ValueError:
                        # ignore non-numeric tokens silently
                        pass

                if parsed_any:
                    started_values = True
                else:
                    # No numeric tokens found:
                    #   - if values haven't started yet, keep scanning
                    #   - if values already started, this is likely the next section -> stop
                    if started_values:
                        break

    except Exception:
        return None

    return vals if vals else None


def auto_find_mesh(sim: Optional[str], root: Path, explicit_mesh: Optional[Path]) -> Path:
    """
    Resolve mesh path:
    - If --mesh is given, use it.
    - Else if --sim is given, use sim_res/<sim>/mesh.mesh under --root.
    - Else autodetect the first *.mesh in CWD.
    """
    if explicit_mesh:
        return explicit_mesh.resolve()

    if sim:
        path = (root / sim / "mesh.mesh")
        if path.exists():
            return path.resolve()
        # try fallback: search under sim folder
        hits = list((root / sim).glob("**/*.mesh"))
        if hits:
            return hits[0].resolve()
        raise SystemExit(f"No mesh found under {root/sim}")

    # bare autodetect in CWD
    hits = list(Path(".").glob("*.mesh"))
    if not hits:
        raise SystemExit("No *.mesh found in current directory. Use --mesh or --sim.")
    if len(hits) > 1:
        print(f"[info] multiple .mesh files found, using: {hits[0]}")
    return hits[0].resolve()


def resolve_interfaces(sim: Optional[str], root: Path,
                       interfaces_arg: Optional[str]) -> Optional[List[float]]:
    """
    Decide interfaces list:
    - If --interfaces is provided: parse comma-separated list or a file path.
    - Else if --sim is provided: try sim_res/<sim>/config.txt section.
    - Else: None.
    """
    if interfaces_arg:
        p = Path(interfaces_arg)
        if p.exists():
            # read numbers from file (one per line or space/comma separated)
            vals: List[float] = []
            txt = p.read_text(encoding="utf-8").replace(",", " ")
            for tok in txt.split():
                try:
                    vals.append(float(tok))
                except ValueError:
                    pass
            return vals if vals else None
        else:
            # comma/space seplist
            vals: List[float] = []
            for tok in interfaces_arg.replace(",", " ").split():
                try:
                    vals.append(float(tok))
                except ValueError:
                    pass
            return vals if vals else None

    if sim:
        cfg = (root / sim / "config.txt")
        vals = parse_interfaces_from_config(cfg)
        return vals

    return None


# --------------------------- plotting utils ---------------------------

def set_axes_equal(ax):
    """Make 3D axes have equal scale for x/y/z."""
    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()
    x_range = abs(x_limits[1] - x_limits[0])
    y_range = abs(y_limits[1] - y_limits[0])
    z_range = abs(z_limits[1] - z_limits[0])
    max_range = max([x_range, y_range, z_range])
    x_middle = np.mean(x_limits)
    y_middle = np.mean(y_limits)
    z_middle = np.mean(z_limits)
    ax.set_xlim3d([x_middle - max_range / 2, x_middle + max_range / 2])
    ax.set_ylim3d([y_middle - max_range / 2, y_middle + max_range / 2])
    ax.set_zlim3d([z_middle - max_range / 2, z_middle + max_range / 2])


def color_wheel(n: int) -> np.ndarray:
    """Return n distinct RGBA colors."""
    cmap = plt.get_cmap("tab20") if n <= 20 else plt.get_cmap("hsv")
    return cmap(np.linspace(0, 1, n))


def build_transition_groups(faces: np.ndarray, doms: np.ndarray) -> Dict[Tuple[int, int], np.ndarray]:
    """
    Group face indices by (front, back) transition.
    Returns dict: (front, back) -> indices (int array into faces).
    """
    groups: Dict[Tuple[int, int], List[int]] = {}
    for idx, (fr, bk) in enumerate(doms):
        key = (int(fr), int(bk))
        groups.setdefault(key, []).append(idx)
    return {k: np.asarray(v, dtype=np.int32) for k, v in groups.items()}


def decimate_indices(idxs: np.ndarray, every: int) -> np.ndarray:
    """Keep every 'every'-th index for speed when plotting huge meshes."""
    if every <= 1:
        return idxs
    return idxs[::every]


# --------------------------- main plotter ---------------------------
def plot_mesh_3d(nodes: np.ndarray,
                 faces: np.ndarray,
                 doms: np.ndarray,
                 interfaces: Optional[List[float]] = None,
                 title: Optional[str] = None,
                 decimate: int = 1,
                 alpha_mesh: float = 0.85,
                 ax=None):
    """
    Render the mesh grouped by (front,back) transitions and optional z-interfaces.
    If 'ax' is provided, draw into that Axes3D. Otherwise create a new figure.
    Returns (fig, ax).
    """
    if ax is None:
        fig = plt.figure(figsize=(9, 7))
        ax = fig.add_subplot(111, projection="3d")
    else:
        fig = ax.figure
        ax.clear()

    # group faces by transition
    groups = build_transition_groups(faces, doms)
    keys = sorted(groups.keys())
    cols = color_wheel(len(keys))

    # build Poly3DCollection with per-face colors
    polys_all = []
    facecolors_all = []
    handles, labels = [], []

    for ci, key in enumerate(keys):
        idxs = decimate_indices(groups[key], decimate)
        if idxs.size == 0:
            continue
        tris = faces[idxs]  # (K,3)
        # build a list of 3D polygons
        polys = [nodes[t] for t in tris]
        polys_all.extend(polys)
        
        # repeat this group's color for its faces
        facecolors_all.extend([cols[ci]] * len(polys))

        # legend handle
        handles.append(Patch(facecolor=cols[ci], edgecolor="k", alpha=alpha_mesh))
        labels.append(f"{key[0]}→{key[1]}")

    if polys_all:
        coll = Poly3DCollection(
            polys_all,
            facecolors=facecolors_all,
            edgecolors="k",
            linewidths=0.15,
            alpha=alpha_mesh,
        )
        ax.add_collection3d(coll)

    # set bounds from nodes
    ax.auto_scale_xyz(nodes[:, 0], nodes[:, 1], nodes[:, 2])
    set_axes_equal(ax)
    # Equal scaling (no squash): use both methods for robustness across Matplotlib versions
    try:
        # Matplotlib >= 3.3
        ax.set_box_aspect((1, 1, 1))   # ensures equal aspect of the 3D box
    except Exception:
        pass

    set_axes_equal(ax)  # keeps data limits spherical

    # draw interface planes
    if interfaces:
        xmn, xmx = ax.get_xlim3d()
        ymn, ymx = ax.get_ylim3d()
        for zi in interfaces:
            Xp, Yp = np.meshgrid([xmn, xmx], [ymn, ymx])
            Zp = np.full_like(Xp, zi, dtype=float)
            plane = ax.plot_surface(Xp, Yp, Zp, rstride=1, cstride=1,
                                    color="0.8", alpha=0.25, linewidth=0, antialiased=False, zorder=3)
            # thin outline
            ax.plot([xmn, xmx, xmx, xmn, xmn], [ymn, ymn, ymx, ymx, ymn], [zi]*5, color="0.4", linewidth=0.6, zorder=4)

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    if title:
        ax.set_title(title)

    # colored legend matching element colors
    if handles:
        leg = ax.legend(handles=handles, labels=labels,
                        title="Domain boundaries (front→back):\n background domain = 0",
                        loc="best", bbox_to_anchor=(1.02, 1.0),
                        borderaxespad=0.0, frameon=True)
        leg.get_frame().set_alpha(0.8)

        # --- centered, multi-line legend title ---
        legend_title = "Domain boundaries (front→back):\nbackground domain = 0"
        leg.set_title(legend_title)
        leg.get_title().set_multialignment("center")  # center each line
        leg.get_title().set_ha("center")              # center the whole block

        # (optional) also center handles/labels block under the title (matplotlib>=3.5-ish)
        try:
            leg._legend_box.align = "center"          # fallback; private API but works widely
        except Exception:
            pass

    return fig, ax


# --------------------------- CLI ---------------------------

def main():
    ap = argparse.ArgumentParser(description="3D plot of HELIOS mesh.mesh with optional interface planes.")
    src = ap.add_mutually_exclusive_group()
    src.add_argument("--sim", type=str, help="Simulation name (uses --root/<sim>/mesh.mesh and config.txt)")
    src.add_argument("--mesh", type=str, help="Path to mesh.mesh")

    ap.add_argument("--root", type=str, default="sim_res", help="Project results root (default: sim_res)")
    ap.add_argument("--interfaces", type=str, default=None,
                    help="Comma-separated z list OR a file containing z values. "
                         "If omitted and --sim is given, attempts to read from config.txt.")
    ap.add_argument("--decimate", type=int, default=1,
                    help="Plot every n-th triangle (speedup for huge meshes).")
    ap.add_argument("--alpha", type=float, default=0.85, help="Mesh face opacity [0..1].")
    ap.add_argument("--title", type=str, default=None, help="Figure title override.")
    ap.add_argument("--save", type=str, default=None, help="Save figure to this path.")
    ap.add_argument("--show", action="store_true", help="Show interactive window.")

    args = ap.parse_args()
    root = Path(args.root)

    mesh_path = auto_find_mesh(args.sim, root, Path(args.mesh) if args.mesh else None)
    interfaces = resolve_interfaces(args.sim, root, args.interfaces)

    nodes, faces, doms = read_mesh_mesh(mesh_path)
    ttl = args.title or (f"{mesh_path.name}" + (f" — {args.sim}" if args.sim else ""))

    fig, _ = plot_mesh_3d(nodes, faces, doms, interfaces=interfaces,
                          title=ttl, decimate=max(1, args.decimate), alpha_mesh=args.alpha)

    if args.save:
        outp = Path(args.save)
        fig.savefig(outp, dpi=200, bbox_inches="tight")
        print(f"[info] saved: {outp}")

    # Show unless running headless
    if args.show or (not args.save):
        plt.show()


if __name__ == "__main__":
    main()
