#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
jobwriter.py — Generate HELIOS job XMLs from a compact config.txt.

This tool reads a user-friendly configuration under sim_data/<SimName>/ and
produces HELIOS job files under sim_res/<SimName>/. It supports isolated and
periodic homogeneous background simulations as well as isolated scatterers in 
layered media backgrounds. It also supports multiple illuminations. Material
permittivities can be provided as literals or via tabulated files.

───────────────────────────────────────────────────────────────────────────────
Project layout & inputs
───────────────────────────────────────────────────────────────────────────────
sim_data/<SimName>/
  └─ config.txt             # textual template describing domain(s), sweep, and incidents
materials/ (optional)
  └─ <material>.txt         # ε(λ) or n,k tables used when referenced by name

Output (default):
sim_res/<SimName>/
  ├─ jobs/
  │   └─ job.1.xml, job.2.xml, ...      # one per wavelength entry
  ├─ Convert.xml                        # mesh converter reference when needed
  └─ angles.txt                         # record of θ,φ for reproducibility

───────────────────────────────────────────────────────────────────────────────
Modes
───────────────────────────────────────────────────────────────────────────────
MODE 0  Isolated homogeneous domain
MODE 1  Periodic homogeneous domain (adds <periodic> block)
MODE 2  Isolated layered medium (expands 'LayeredM#' tokens into <domainLayers>)

───────────────────────────────────────────────────────────────────────────────
Materials & permittivity
───────────────────────────────────────────────────────────────────────────────
- Literal ε: tokens like "2.25+0.0j" are treated as constant complex εr.
- Tabulated files: DATA_DIR/<name>.txt is read when a material token names a file.
  Two formats are accepted:
    - [λ[nm], Re(ε), Im(ε)]
    - [λ[nm], n, k]    (enable via env MATERIALS_COLUMNS=nk)
  Linear interpolation is used across λ; endpoints are linearly extrapolated.

Environment variables:
  DATA_DIR              Override default materials directory (../materials)
  MATERIALS_COLUMNS     "nk" to interpret tables as (n,k) instead of (Re(ε),Im(ε))

───────────────────────────────────────────────────────────────────────────────
Incidence & sweeps
───────────────────────────────────────────────────────────────────────────────
- Angles: θ,φ may be single values or ranges "start end step" (inclusive if exact).
- Polarization: s/p (TE/TM) or circular (LCP/RCP); complex polarization vectors
  are emitted when required by the chosen basis.
- Each θ,φ pair yields a <job> (with optional <subjob>s for multiple illuminations).

───────────────────────────────────────────────────────────────────────────────
CLI synopsis (project-native form)
───────────────────────────────────────────────────────────────────────────────
jobwriter.py MODE COMSOL SimName
  - Reads sim_data/<SimName>/config.txt
  - Writes to sim_res/<SimName>/{jobs/, angles.txt}

───────────────────────────────────────────────────────────────────────────────
Notes
───────────────────────────────────────────────────────────────────────────────
- In layered mode, two extra lines in config.txt enumerate layer materials and
  their interface z-positions per stack. Occurrences of 'LayeredM#' within the
  domain line are expanded in order to the corresponding stacks.
- For periodic homogeneous jobs, the <periodic> block carries in-plane lattice
  vectors and the Bloch wavevector consistent with the incidence (λ and ε).

"""

from __future__ import division
import os
import sys
from pathlib import Path
from typing import List, Tuple, Dict, Union
import numpy as np
from numpy import pi, sin, cos, cross, array, sqrt, real, imag
from scipy import interpolate
import re

# -----------------------------------------------------------------------------
# Paths / constants
# -----------------------------------------------------------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR   = SCRIPT_DIR.parent

# New standard locations
SIM_DATA_ROOT = ROOT_DIR / 'sim_data'    # input configs+meshes live here
SIM_RES_ROOT  = ROOT_DIR / 'sim_res'     # outputs go here

# MATERIALS_DIR overrides default '../materials' for ε(λ) tables (handy on clusters/shared volumes).
DATA_DIR = Path(os.getenv('MATERIALS_DIR', ROOT_DIR / 'materials')).expanduser().resolve()

EPS = 1e-12  # numerical zero

# -----------------------------------------------------------------------------
# Small helpers
# -----------------------------------------------------------------------------
def _parse_eps_literal(s: str):
    """
    Try to parse a permittivity literal from string `s`.
    Accepted forms (whitespace allowed):
      - Real:            1      |  2.25  |  -3.1e-2
      - Complex (rect):  2.25+0.1i   |  2.25-0.1i   |  2.25+0.1j
      - Complex (tuple): (2.25,0.1)  |  [2.25, 0.1] | {2.25, 0.1}

    Returns complex(...) on success, or None if `s` is not a literal.
    """
    t = s.strip()

    # 1) Tuple-like: (a,b) / [a,b] / {a,b}
    m = re.match(
        r'^[\(\[\{]\s*([+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?)\s*[,;]\s*'
        r'([+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?)\s*[\)\]\}]$',
        t
    )
    if m:
        return complex(float(m.group(1)), float(m.group(2)))

    # 2) a+bi / a-bi (allow 'i' or 'j')
    tj = t.replace('I', 'i').replace('i', 'j')
    # Only allow characters that Python's complex() understands
    if re.fullmatch(r'[ \t\+\-\d\.\(\)eEjJ]+', tj):
        try:
            return complex(tj)
        except ValueError:
            pass

    # 3) Pure real
    try:
        return complex(float(t), 0.0)
    except ValueError:
        return None

def lin_or_range(s: str) -> np.ndarray:
    """
    Parse one number or 'start stop step' -> numpy array (inclusive end if exact).
    """
    parts = s.split()
    if len(parts) == 1:
        return np.array([float(parts[0])], dtype=float)
    if len(parts) == 3:
        a, b, c = map(float, parts)
        if c == 0:
            raise ValueError("step must be nonzero")
        n = round((b - a) / c)
        arr = a + c * np.arange(int(n) + 1, dtype=float)
        if arr.size == 0 or abs(arr[-1] - b) > 1e-10 * max(1.0, abs(b)):
            # ensure inclusion if close
            if (b - a) / c > 0:
                arr = a + c * np.arange(int((b - a) / c) + 1 + 1, dtype=float)
                arr = arr[arr <= b + 1e-12]
        return arr
    raise ValueError("Expected 1 number or 'start stop step'")

def clean_lines(path: Path) -> List[str]:
    """
    Load config, strip blank lines and comments ('#' at line start).
    Keep section order identical to user templates.
    """
    out = []
    with open(path, 'r', encoding='utf-8-sig', errors='ignore') as f:
        for line in f:
            t = line.strip()
            if not t:
                continue
            if t.startswith('#'):
                continue
            out.append(t)
    return out

def split_layered_list(line: str) -> List[List[str]]:
    """
    Parse layered materials line: stacks are comma-separated; within a stack,
    materials are space-separated; return list of stacks (each a list of tokens).
    Example: "1 SiO2 1, Ag Au" -> [["1","SiO2","1"], ["Ag","Au"]]
    """
    stacks = []
    for seg in line.split(','):
        tokens = [tok for tok in seg.strip().split() if tok]
        if tokens:
            stacks.append(tokens)
    return stacks

def split_interfaces_list(line: str) -> List[List[float]]:
    """
    Parse layered interfaces line: stacks are comma-separated; within a stack,
    numbers are space-separated; return list of stacks (each a list of floats).
    Example: "0 150, 20" -> [[0.0, 150.0], [20.0]]
    """
    stacks = []
    for seg in line.split(','):
        toks = [tok for tok in seg.strip().split() if tok]
        if toks:
            stacks.append([float(x) for x in toks])
    return stacks

def get_eps(material: str, lam_nm: np.ndarray, use_spline: bool = False) -> np.ndarray:
    """
    Return ε(λ) for a material at wavelengths lam_nm [nm].

    Accepted file formats at DATA_DIR/<material>.txt (any number of header lines ok):
      - λ  Re(ε)  Im(ε)           (default)
      - λ  n      k               (set env MATERIALS_COLUMNS=nk to convert to ε)

    The parser ignores blank lines, lines starting with '#', and any non-numeric tokens.
    Commas are treated as whitespace, so CSV-like lines are fine.
    """
    val = _parse_eps_literal(material)
    if val is not None:
        # constant ε over wavelength
        return np.full_like(lam_nm, val, dtype=complex)


    path = DATA_DIR / f"{material}.txt"
    if not path.exists():
        raise FileNotFoundError(f"Material file not found: {path}")

    rows = []
    with open(path, 'r', encoding='utf-8-sig', errors='ignore') as f:
        for line in f:
            # strip comments / blanks
            s = line.strip()
            if not s or s.startswith('#'):
                continue
            # treat commas like spaces
            s = s.replace(',', ' ')
            # keep only numeric tokens
            nums = []
            for tok in s.split():
                try:
                    nums.append(float(tok))
                except ValueError:
                    # ignore non-numeric tokens in header-like lines
                    continue
            if len(nums) >= 3:
                # use first three numeric tokens: λ, col2, col3
                rows.append(nums[:3])

    if not rows:
        raise RuntimeError(f"No numeric data found in {path}")

    arr = np.asarray(rows, dtype=float)
    wav  = arr[:, 0]
    col2 = arr[:, 1]
    col3 = arr[:, 2]

    # Decide whether columns are eps_r/eps_i or n/k
    mode = os.getenv('MATERIALS_COLUMNS', 'eps').lower()  # 'eps' or 'nk'
    if mode == 'nk':
        n = col2
        k = col3
        eps = (n + 1j * k) ** 2
    else:
        eps = col2 + 1j * col3

    # Interpolate ε(λ)
    # - Default: linear (previous behavior)
    # - With --spline: cubic spline, similar to MATLAB spline(...)
    if use_spline:
        rfun = interpolate.CubicSpline(wav, np.real(eps), extrapolate=True)
        ifun = interpolate.CubicSpline(wav, np.imag(eps), extrapolate=True)
    else:
        rfun = interpolate.interp1d(
            wav, np.real(eps), kind='linear',
            bounds_error=False, fill_value='extrapolate'
        )
        ifun = interpolate.interp1d(
            wav, np.imag(eps), kind='linear',
            bounds_error=False, fill_value='extrapolate'
        )
    return rfun(lam_nm) + 1j * ifun(lam_nm)

def build_incident(theta_deg: float, phi_deg: float, pol: str):
    """
    Return ('planewave', kdir, polvec) supporting pol in {s,p,TE,TM,LCP,RCP}.
    """
    pol = pol.lower()
    theta = theta_deg * pi / 180.0
    phi   = phi_deg   * pi / 180.0
    kdir = array([sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)], dtype=float)
    kdir[abs(kdir) < EPS] = 0.0
    ex = array([sin(phi), -cos(phi), 0.0], dtype=float)  # s/TE
    ey = cross(kdir, ex)                                 # p/TM
    if pol in ('s','te'):
        polvec = ex
    elif pol in ('p','tm'):
        polvec = ey
    elif pol == 'lcp':
        polvec = (ex + 1j*ey) / np.sqrt(2.0)
    elif pol == 'rcp':
        polvec = (ex - 1j*ey) / np.sqrt(2.0)
    else:
        raise ValueError("Polarization must be one of: s, p, TE, TM, LCP, RCP")
    polvec[abs(polvec) < EPS] = 0.0
    return ('planewave', kdir, polvec)

def _unit_vec_from_angles(theta_deg: float, phi_deg: float) -> np.ndarray:
    th = theta_deg * pi / 180.0
    ph = phi_deg   * pi / 180.0
    v = np.array([sin(th)*cos(ph), sin(th)*sin(ph), cos(th)], dtype=float)
    v[abs(v) < EPS] = 0.0
    return v

def build_dipole(loc_xyz: Tuple[float,float,float],
                 pol: Union[Tuple[float,float,float], Tuple[float,float]],
                 kind: str = "vec"):
    """
    Return ('dipole', location[3], polarization[3]).
    - kind='vec'  => pol = (px,py,pz)
    - kind='ang'  => pol = (theta_deg, phi_deg)  (turns into unit vector)
    """
    if kind == "ang":
        pv = _unit_vec_from_angles(pol[0], pol[1])
    else:
        pv = np.array(pol, dtype=float)
    pv[abs(pv) < EPS] = 0.0
    return ('dipole', np.array(loc_xyz, dtype=float), pv)

def emit_dipole_block(dip: tuple) -> str:
    kind, r0, pv = dip
    if kind != 'dipole':
        raise ValueError("emit_dipole_block expects a 'dipole' tuple")
    s  = "\t<dipole>\n"
    s += f"\t\t<location> {r0[0]:.12f} {r0[1]:.12f} {r0[2]:.12f} </location>\n"
    s += f"\t\t<polarization> {pv[0]:.12f} {pv[1]:.12f} {pv[2]:.12f} </polarization>\n"
    s += "\t</dipole>\n"
    return s

def _max_sub_index(sub_map: Dict[int, List[Tuple]]) -> int:
    return max(sub_map.keys(), default=1)

def _apply_to_targets(subs: List[int], incs: List[Tuple],
                      main_incs: List[Tuple],
                      sub_map: Dict[int, List[Tuple]],
                      mode: str) -> None:
    """
    Apply 'incs' to each target in 'subs' according to mode:
      - mode == 'overwrite': replace targeted subjob contents.
      - mode == 'add': append to existing.
    - 'subs' may contain 1 ('main') or subjob indices >= 2.
    - if target index > last existing, append at (last+1) instead
    """
    for target in subs:
        if target == 1:
            if mode == 'overwrite':
                main_incs.clear()
            main_incs.extend(incs)
        else:
            last = _max_sub_index(sub_map)
            if mode == 'overwrite':
                if target <= last and target in sub_map:
                    sub_map[target] = list(incs)
                else:
                    # Target is beyond current range -> append at (last+1)
                    new_idx = last + 1
                    sub_map[new_idx] = list(incs)
            else:  # mode == 'add'
                if target <= last and target in sub_map:
                    sub_map[target].extend(incs)
                else:
                    # Target is beyond current range -> append at (last+1)
                    new_idx = last + 1
                    sub_map.setdefault(new_idx, []).extend(incs)

# -----------------------------------------------------------------------------
# XML emitters
# -----------------------------------------------------------------------------
def emit_convert(mesh_name: str = 'mesh.mphtxt') -> str:
    """
    Convert.xml that references the local mesh file inside the same output folder.
    """
    base = '.'.join(mesh_name.split('.')[:-1]) or mesh_name
    return (
        "<conversion>\n"
        f"\t<mesh>\n\t{mesh_name}\n\t</mesh>\n\n"
        f"\t<out>\n\t{base}\n\t</out>\n"
        "</conversion>\n"
    )

def emit_incident_block(inc: tuple, inci_domain: int) -> str:
    """
    <planewave> ... </planewave> from ('planewave', kdir, polvec).
    Supports both real (s/p/TE/TM) and complex (LCP/RCP) polarization vectors.
    """
    kind, kdir, polvec = inc
    if kind != 'planewave':
        raise ValueError("emit_incident_block expects a 'planewave'")

    s = "\t<planewave>\n"
    s += f"\t\t<domain> {inci_domain} </domain>\n"
    s += f"\t\t<propagation> {kdir[0]:.12f} {kdir[1]:.12f} {kdir[2]:.12f} </propagation>\n"

    pv = np.asarray(polvec)
    if np.iscomplexobj(pv):
        s += (
            "\t\t<polarization> "
            f"({pv[0].real:.12f},{pv[0].imag:.12f}) "
            f"({pv[1].real:.12f},{pv[1].imag:.12f}) "
            f"({pv[2].real:.12f},{pv[2].imag:.12f}) </polarization>\n"
        )
    else:
        s += f"\t\t<polarization> {pv[0]:.12f} {pv[1]:.12f} {pv[2]:.12f} </polarization>\n"

    s += "\t</planewave>\n"
    return s

def job_iso_hom(job_idx: int, mesh_name: str, wl: float,
                eps_domains: np.ndarray,
                main_incs: List[Tuple], subjob_map: Dict[int, List[Tuple]],
                inci_domain: int) -> str:
    """
    <job> for isolated homogeneous.
    """
    s = "<job>\n"
    s += f"\t<label> out/sim{job_idx}_inc_1 </label>\n"
    s += f"\t<mesh> {mesh_name} </mesh>\n"
    s += f"\t<wavelength> {wl:.12f} </wavelength>\n"
    for eps in eps_domains:
        s += f"\t<domain> ({eps.real:.12f},{eps.imag:.12f}) (1.0,0.0) </domain>\n"
    # Emit all main-job illuminations in order
    for inc in main_incs:
        if inc[0] == 'planewave':
            s += emit_incident_block(inc, inci_domain)
        else:
            s += emit_dipole_block(inc)
    s += "</job>\n"
    # Explicit subjobs (keys are the desired *_inc_<N>)
    for sub_idx in sorted(subjob_map.keys()):
        s += "\n<subjob>\n"
        s += f"\t<label> out/sim{job_idx}_inc_{sub_idx} </label>\n"
        for inc in subjob_map[sub_idx]:
            if inc[0] == 'planewave':
                s += emit_incident_block(inc, inci_domain)
            else:
                s += emit_dipole_block(inc)
        s += "</subjob>\n"
    return s

def job_per_hom(job_idx: int, mesh_name: str, wl: float,
                eps_domains: np.ndarray,
                main_incs: List[Tuple], subjob_map: Dict[int, List[Tuple]],
                inci_domain: int,
                lattice2: List[List[float]]) -> str:
    """
    <job> for periodic homogeneous: add <periodic> with lattice and Bloch.
    """
    # Bloch from first incident direction
    kdir = main_incs[0][1]
    k_mag = 2.0 * pi * sqrt(max(0.0, real(eps_domains[inci_domain]))) / wl
    bloch = k_mag * kdir
    s = "<job>\n"
    s += f"\t<label> out/sim{job_idx}_inc_1 </label>\n"
    s += "\t<periodic>\n"
    s += f"\t\t<lattice> {lattice2[0][0]} {lattice2[0][1]} {lattice2[0][2]} </lattice>\n"
    s += f"\t\t<lattice> {lattice2[1][0]} {lattice2[1][1]} {lattice2[1][2]} </lattice>\n"
    s += f"\t\t<bloch> {bloch[0]:.12f} {bloch[1]:.12f} {bloch[2]:.12f} </bloch>\n"
    s += "\t</periodic>\n"
    s += f"\t<mesh> {mesh_name} </mesh>\n"
    s += f"\t<wavelength> {wl:.12f} </wavelength>\n"
    for eps in eps_domains:
        s += f"\t<domain> ({eps.real:.12f},{eps.imag:.12f}) (1.0,0.0) </domain>\n"
    for inc in main_incs:
        if inc[0] == 'planewave':
            s += emit_incident_block(inc, inci_domain)
        else:
            s += emit_dipole_block(inc)
    s += "</job>\n"
    for sub_idx in sorted(subjob_map.keys()):
        s += "\n<subjob>\n"
        s += f"\t<label> out/sim{job_idx}_inc_{sub_idx} </label>\n"
        for inc in subjob_map[sub_idx]:
            if inc[0] == 'planewave':
                s += emit_incident_block(inc, inci_domain)
            else:
                s += emit_dipole_block(inc)
        s += "</subjob>\n"
    return s

def job_iso_layered(job_idx: int, mesh_name: str, wl: float,
                    domain_tokens: List[str],
                    eps_hom_this_wl: Dict[str, complex],
                    stacks_eps_this_wl: List[List[complex]],
                    stacks_interfaces: List[List[float]],
                    main_incs: List[Tuple], 
                    subjob_map: Dict[int, List[Tuple]],
                    inci_domain: int) -> str:
    """
    <job> for isolated layered:
      - Original domain token order.
      - For 'LayeredM#' emit a <domainLayers> with layers & (optional) interfaces.
      - For plain materials (e.g., "1", "Ag"), emit a <domain>.
    """
    s = "<job>\n"
    s += f"\t<label> out/sim{job_idx}_inc_1 </label>\n"
    s += f"\t<mesh> {mesh_name} </mesh>\n"
    s += f"\t<wavelength> {wl:.12f} </wavelength>\n"

    # Map LayeredM index to stack position (1-based in token -> 0-based in list)
    lay_counter = 0
    for tok in domain_tokens:
        if tok.startswith("LayeredM"):
            # which layered stack?
            lay_counter += 1
            stack_idx = lay_counter - 1
            layers = stacks_eps_this_wl[stack_idx] if stack_idx < len(stacks_eps_this_wl) else []
            zifs   = stacks_interfaces[stack_idx]  if stack_idx < len(stacks_interfaces)  else []
            s += "\t<domainLayers>\n"
            for li, eps in enumerate(layers):
                s += f"\t\t<layer> ({eps.real:.12f},{eps.imag:.12f}) (1.0,0.0) </layer>\n"
                if li < len(layers) - 1 and li < len(zifs):
                    s += f"\t\t<interface> {zifs[li]} </interface>\n"
            s += "\t</domainLayers>\n"
        else:
            eps = eps_hom_this_wl[tok]
            s += f"\t<domain> ({eps.real:.12f},{eps.imag:.12f}) (1.0,0.0) </domain>\n"

    for inc in main_incs:
        if inc[0] == 'planewave':
            s += emit_incident_block(inc, inci_domain)
        else:
            s += emit_dipole_block(inc)
    s += "</job>\n"
    for sub_idx in sorted(subjob_map.keys()):
        s += "\n<subjob>\n"
        s += f"\t<label> out/sim{job_idx}_inc_{sub_idx} </label>\n"
        for inc in subjob_map[sub_idx]:
            if inc[0] == 'planewave':
                s += emit_incident_block(inc, inci_domain)
            else:
                s += emit_dipole_block(inc)
        s += "</subjob>\n"
    return s

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
def main():
    # Optional CLI flag:
    #   --spline  -> use cubic spline interpolation for ε(λ) instead of linear
    use_spline = '--spline' in sys.argv[1:]
    argv = [a for a in sys.argv if a != '--spline']

    # Support both:
    #  - New short form: jobwriter.py MODE COMSOL SimName [--spline]
    #  - Legacy form   : jobwriter.py MODE COMSOL CONFIG SimName [SimName2 ...] [--spline]
    if len(argv) < 4:
        print(f"Usage (new):   {argv[0]} MODE(0/1/2) COMSOL(0/1) SimName [--spline]")
        print(f"Usage (legacy):{argv[0]} MODE(0/1/2) COMSOL(0/1) CONFIG path_folder1 [path_folder2 ...] [--spline]")
        print("  MODE 0: isolated homogeneous")
        print("  MODE 1: periodic  homogeneous")
        print("  MODE 2: isolated layered")
        sys.exit(1)

    MODE = int(argv[1])

    # Decide which calling convention we are in
    if len(argv) == 4:
        # --- NEW SHORT FORM ---
        sim_name = argv[3]
        input_dir  = (SIM_DATA_ROOT / sim_name).resolve()
        config_path = input_dir / 'config.txt'
        folders = [sim_name]
        base_out = (SIM_RES_ROOT / sim_name).resolve()
    else:
        # --- LEGACY FORM ---
        config_path = Path(argv[3]).resolve()
        folders = argv[4:]
        base_out = (SIM_RES_ROOT / sim_name).resolve()

    # Load and sanitize config lines
    lines = clean_lines(config_path)
    it = iter(lines)

    # 1) Wavelengths
    wl_line = next(it)
    wavelengths = lin_or_range(wl_line)

    # 2) Domain tokens
    domain_line = next(it)
    domain_tokens = domain_line.split()

    # 3) MODE-specific layered blocks
    layered_mat_stacks = []
    layered_if_stacks  = []
    if MODE == 2:
        mats_line = next(it)  # e.g. "1 SiO2 1, Ag Au"
        layered_mat_stacks = split_layered_list(mats_line)
        if_line  = next(it)  # e.g. "0 150, 20"
        layered_if_stacks  = split_interfaces_list(if_line)
        while len(layered_if_stacks) < len(layered_mat_stacks):
            layered_if_stacks.append([])

    # 4) Incident domain
    inci_domain = int(next(it))

    # 5) Polarization
    pol = next(it).strip()

    # 6) Theta / 7) Phi
    thetas = lin_or_range(next(it))
    phis   = lin_or_range(next(it))

    # 8) Periodic params if MODE==1
    lattice2 = None
    if MODE == 1:
        px, py, _, _, _ = map(float, next(it).split())
        lattice2 = [[px, 0.0, 0.0], [0.0, py, 0.0]]

    # 9) Optional per-subjob illumination directives
    def _parse_sub_ranges(spec: str) -> List[int]:
        # spec like "2-8,10,13,15-20" or "main"
        out: List[int] = []
        for chunk in spec.split(','):
            c = chunk.strip().lower()
            if c == 'main':
                out.append(1)  # main job index
                continue
            if '-' in c:
                a, b = c.split('-', 1)
                out.extend(range(int(a), int(b)+1))
            else:
                out.append(int(c))
        return sorted(set(out))

    def _parse_sweep_kv(tokens: List[str]) -> Dict[str, str]:
        """
        Parse tokens like ["x=a","b","c","y=d","e","f","theta=10","20","5","subs=2-3,10"]
        into {"x":"a b c", "y":"d e f", "theta":"10 20 5", "subs":"2-3,10"}.
        Stops at '#' (comment).
        """
        out: Dict[str, str] = {}
        key = None
        acc: List[str] = []
        for tok in tokens:
            if tok.startswith('#'):
                break
            if '=' in tok:
                # flush previous
                if key is not None:
                    out[key] = ' '.join(acc).strip()
                key, val = tok.split('=', 1)
                key = key.strip().lower()
                acc = [val.strip()]
            else:
                if key is not None and tok.strip():
                    acc.append(tok.strip())
                # else: stray token before first key, ignore
        if key is not None:
            out[key] = ' '.join(acc).strip()
        return out

    # Collect user-directed incidents
    # - main_incs: list of illuminations for <job> (index==1)
    # - sub_map  : dict subjob_index -> list of illuminations
    main_incs: List[Tuple] = []
    sub_map : Dict[int, List[Tuple]] = {}
    remaining = list(it)  # whatever is left in the config
    # Illumination combination mode: 'overwrite' or 'add' (default)
    illum_mode = 'add'
    # Collect directives first; we will apply them after creating the sweep baseline.
    # Each directive: (subs: List[int], incs: List[Tuple]])
    pending_directives: List[Tuple[List[int], List[Tuple]]] = []
    for line in remaining:
        t = line.strip()
        if not t:
            continue
        # Global combination rule for sweep vs. directives
        #   INC-MODE add        -> keep sweep AND directives (default)
        #   INC-MODE overwrite  -> directives replace sweep
        if t.upper().startswith('INC-MODE'):
            _, mode = t.split(None, 1)
            mode = mode.strip().lower()
            if mode not in ('add', 'overwrite'):
                raise SystemExit("INC-MODE must be 'add' or 'overwrite'")
            illum_mode = mode
            continue
        # Syntax: PW theta phi pol subs=2-8,10 or subs=main
        if t.upper().startswith('PW '):
            _, th, ph, ppol, *rest = t.split()
            subs = [1]  # default to main if no subs= provided
            for tok in rest:
                if tok.lower().startswith('subs='):
                    subs = _parse_sub_ranges(tok.split('=',1)[1])
            inc = build_incident(float(th), float(ph), ppol)
            pending_directives.append((subs, [inc]))
            continue
        #   DIP x y z polang=theta,phi  subs=...   (angles -> unit vector)
        #   DIP x y z polvec=px,py,pz   subs=...   (direct vector)
        if t.upper().startswith('DIP '):
            parts = t.split()
            _, x, y, z, *rest = parts
            polmode = None
            polval  = None
            subs = []
            for tok in rest:
                tl = tok.lower()
                if tl.startswith('polang='):
                    a,b = tl.split('=',1)[1].split(',',1)
                    polmode = 'ang'; polval = (float(a), float(b))
                elif tl.startswith('polvec='):
                    a,b,c = tl.split('=',1)[1].split(',',2)
                    polmode = 'vec'; polval = (float(a), float(b), float(c))
                elif tl.startswith('subs='):
                    subs = _parse_sub_ranges(tl.split('=',1)[1])
            if polmode is None or not subs:
                raise SystemExit("DIP requires polang=θ,φ OR polvec=px,py,pz AND subs=...")
            dip = build_dipole((float(x), float(y), float(z)), polval, polmode)
            pending_directives.append((subs, [dip]))
            continue
        #   DIP-SWEEP x=a b c y=d e f z=g h i theta=j k l phi=m n o subs=2-8,10
        if t.upper().startswith('DIP-SWEEP'):
            # Robustly parse key=value groups where value can contain spaces (e.g., "x=a b c")
            toks = t.split()[1:]
            kv = _parse_sweep_kv(toks)
            # Validate required keys
            required = {'x','y','z','theta','phi','subs'}
            missing = [k for k in required if k not in kv]
            if missing:
                raise SystemExit(f"DIP-SWEEP missing keys: {', '.join(missing)}")
            xs = lin_or_range(kv['x'])
            ys = lin_or_range(kv['y'])
            zs = lin_or_range(kv['z'])
            ths = lin_or_range(kv['theta'])
            phs = lin_or_range(kv['phi'])
            subs = _parse_sub_ranges(kv['subs'])
            # Grouped set: every generated dipole goes into each target subjob
            dips: List[Tuple] = []
            for X in xs:
                for Y in ys:
                    for Z in zs:
                        for TH in ths:
                            for PH in phs:
                                dips.append(build_dipole((float(X),float(Y),float(Z)), (float(TH),float(PH)), 'ang'))
            pending_directives.append((subs, dips))
            continue

    # Collect unique plain (homogeneous) materials from domain_tokens
    plain_materials = [tok for tok in domain_tokens if not tok.startswith("LayeredM")]
    uniq_plain = []
    for m in plain_materials:
        if m not in uniq_plain:
            uniq_plain.append(m)

    # For layered mode, also collect materials used inside layered stacks
    layered_materials = []
    if MODE == 2:
        for stack in layered_mat_stacks:
            for m in stack:
                if m not in layered_materials:
                    layered_materials.append(m)
    # Merge for table building (but keep maps separate)
    all_materials = list(dict.fromkeys(uniq_plain + layered_materials))  # preserve order & unique

    # Interpolate ε(λ) for each material name
    eps_table: Dict[str, np.ndarray] = {}
    for name in all_materials:
        eps_table[name] = get_eps(name, wavelengths, use_spline=use_spline)

    # ------ Prepare default planewave set from theta/phi/pol (back-compatibility) ------
    # Behavior: If user provided any PW/DIP directives, we auto-spawn subjobs
    # from the sweep, then we only place (overwrite or add) the provided any PW/DIP 
    # illuminations in the corresponding subjobs.
    
    # Build baseline subjobs from θ–φ sweep FIRST
    sweep_incs = [build_incident(float(t), float(p), pol) for p in phis for t in thetas]
    if not sweep_incs:
        raise SystemExit("No incidence angles provided.")
    # main job gets the first; subjobs 2..N get the rest
    main_incs.append(sweep_incs[0])
    for k, inc in enumerate(sweep_incs[1:], start=2):
        sub_map.setdefault(k, []).append(inc)

    # THEN apply all directives with the chosen mode
    for subs, incs in pending_directives:
        _apply_to_targets(subs, incs, main_incs, sub_map, illum_mode)

    # Emit per folder
    for folder in folders:
        base = base_out
        jobs = base / 'jobs'
        jobs.mkdir(parents=True, exist_ok=True)

        # Convert.xml (local mesh)
        with open(base / 'Convert.xml', 'w') as f:
            f.write(emit_convert('mesh.mphtxt'))

        # angles.txt
        with open(base / 'angles.txt', 'w') as f:
            for phi in phis:
                for theta in thetas:
                    f.write(f"{theta:g}\t{phi:g}\t{pol}\n")

        # Write job.N.xml with mesh_name='mesh.mphtxt'
        for i, wl in enumerate(wavelengths, start=1):
            wl = float(wl)
            if MODE in (0, 1):
                eps_domains = np.array([eps_table[m][i-1] for m in domain_tokens], dtype=complex)
                if MODE == 0:
                    job = job_iso_hom(i, 'mesh.mesh', wl, eps_domains, main_incs, sub_map, inci_domain)
                else:
                    job = job_per_hom(i, 'mesh.mesh', wl, eps_domains, main_incs, sub_map, inci_domain, lattice2)
            else:
                eps_hom_this_wl = {m: eps_table[m][i-1] for m in domain_tokens if not m.startswith('LayeredM')}
                stacks_eps_this_wl = [[eps_table[m][i-1] for m in stack] for stack in layered_mat_stacks]
                job = job_iso_layered(
                    i, 'mesh.mesh', wl,
                    domain_tokens,
                    eps_hom_this_wl,
                    stacks_eps_this_wl,
                    layered_if_stacks,
                    main_incs, sub_map, inci_domain
                )
            with open(jobs / f"job.{i}.xml", 'w') as f:
                f.write(job)

        print(f"[{folder}] Wrote {len(wavelengths)} job files + angles.txt + Convert.xml to {base}")

if __name__ == '__main__':
    main()
