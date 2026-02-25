#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
meshconvert.py — Convert COMSOL *.mphtxt or Gmsh *.msh (v4 ASCII) meshes to 
                 HELIOS 'mesh.mesh' format.

This converter reads a text mesh exported from COMSOL or Gmsh, reconstructs 
surface triangles, assigns (front, back) domain IDs from neighboring 
tetrahedra, and writes the compact mesh.mesh file expected by HELIOS. Optional 
periodic cropping limits faces to a primitive cell and excludes a cut region.

───────────────────────────────────────────────────────────────────────────────
Inputs & outputs
───────────────────────────────────────────────────────────────────────────────
Input:
  - COMSOL text mesh (*.mphtxt) with sections:
      'Mesh vertex coordinates', '# number of elements', etc.
  - Gmsh MSH v4 ASCII (*.msh): parsed using $Nodes / $Elements / $Entities blocks.
  - Optional periodic parameters (px, py, cx, cy, zcut) to crop to a window
    centered at (cx,cy) with periods (px,py) and to exclude z < zcut.
NOTE: COMSOL and Gmsh inputs support 3D tetrahedral meshes only.

Output:
  mesh.mesh
    Nnodes
    x y z                (per node)
    Nfaces
    i j k  front back    (1-based indices; canonicalized orientation)

───────────────────────────────────────────────────────────────────────────────
Periodic window (optional)
───────────────────────────────────────────────────────────────────────────────
- A triangle is kept only if its centroid lies strictly inside the in-plane
  window and z > zcut. This avoids duplicated faces on periodic boundaries and
  removes unwanted substrate/bulk regions.
- Periodic parameters can be supplied explicitly on the CLI or auto-detected
  from sim_data/<SimName>/config.txt when using the project short form.

───────────────────────────────────────────────────────────────────────────────
Workflow
───────────────────────────────────────────────────────────────────────────────
1) Parse *.mphtxt or *.msh: load nodes, surface triangles and tetrahedra and 
   tags.
2) Build a face->adjacent-tet map; deduce (front, back) domain per triangle by
   comparing a face normal with vectors to neighboring tet centers.
3) Optionally crop by (px,py,cx,cy,zcut).
4) Write mesh.mesh with canonical (front <= back) orientation.

───────────────────────────────────────────────────────────────────────────────
CLI synopsis (selected)
───────────────────────────────────────────────────────────────────────────────
Project short form:
  meshconvert.py SimName
    - Finds sim_data/<SimName>/*.mphtxt (single file expected)
    - Optionally reads (px,py,cx,cy,zcut) from config.txt
    - Writes sim_res/<SimName>/mesh.mesh (and copies source for provenance)

Legacy forms (still supported):
  meshconvert.py px py cx cy zcut
  meshconvert.py file.mphtxt px py cx cy zcut
  meshconvert.py file.mphtxt           # isolated (non-periodic)

───────────────────────────────────────────────────────────────────────────────
Notes & assumptions
───────────────────────────────────────────────────────────────────────────────
- COMSOL text meshes are supported; section markers are used for robust parsing.
- Gmsh MSH v4 ASCII is also supported.
- Domain IDs are taken from tetrahedral region tags; faces shared by two regions
  get (front, back) according to outward normal orientation.
- The writing order (nodes -> faces) matches HELIOS reader expectations exactly.

"""

from __future__ import division
import os
import sys
import shutil
from pathlib import Path
from typing import List, Tuple, Dict, Optional

# --------------------------------------------------------------------------------------
# Project-aware paths
# --------------------------------------------------------------------------------------

SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR   = SCRIPT_DIR.parent
SIM_DATA_ROOT = ROOT_DIR / 'sim_data'   # inputs (config.txt + *.mphtxt)
SIM_RES_ROOT  = ROOT_DIR / 'sim_res'    # outputs (mesh.mphtxt + mesh.mesh)

# --------------------------------------------------------------------------------------
# Small vector utilities
# --------------------------------------------------------------------------------------

def tri_center(a, b, c):
    """Small geometry helper."""
    return [(a[0] + b[0] + c[0]) / 3.0,
            (a[1] + b[1] + c[1]) / 3.0,
            (a[2] + b[2] + c[2]) / 3.0]

def tet_center(a, b, c, d):
    """Small geometry helper."""
    return [(a[0] + b[0] + c[0] + d[0]) / 4.0,
            (a[1] + b[1] + c[1] + d[1]) / 4.0,
            (a[2] + b[2] + c[2] + d[2]) / 4.0]

def cross(u, v):
    """Small geometry helper."""
    return [u[1]*v[2] - u[2]*v[1],
            u[2]*v[0] - u[0]*v[2],
            u[0]*v[1] - u[1]*v[0]]

def dot(u, v):
    """Small geometry helper."""
    return u[0]*v[0] + u[1]*v[1] + u[2]*v[2]

def face_normal(a, b, c):
    """Small geometry helper."""
    ab = [b[0]-a[0], b[1]-a[1], b[2]-a[2]]
    bc = [c[0]-b[0], c[1]-b[1], c[2]-b[2]]
    return cross(ab, bc)

def test_edge_periodic(nodes, tri, PP, eps=0.1):
    """Small geometry helper."""
    a, b, c = nodes[tri[0]], nodes[tri[1]], nodes[tri[2]]
    FC = tri_center(a, b, c)
    px, py, cx, cy, zcut = PP
    if not (cx - px/2.0 + eps < FC[0] < cx + px/2.0 - eps):
        return False
    if not (cy - py/2.0 + eps < FC[1] < cy + py/2.0 - eps):
        return False
    if zcut >= FC[2] - eps:
        return False
    return True

# --------------------------------------------------------------------------------------
# Text parsing helpers for .mphtxt
# --------------------------------------------------------------------------------------

def first_int_in_line(s: str):
    """Extract first/sequence of integers/floats from a whitespace-separated line."""
    for tok in s.strip().split():
        try:
            return int(tok)
        except ValueError:
            continue
    return None

def ints_from_line(s: str) -> List[int]:
    """Extract first/sequence of integers/floats from a whitespace-separated line."""
    return [int(x) for x in s.split()]

def floats_from_line(s: str) -> List[float]:
    """Extract first/sequence of integers/floats from a whitespace-separated line."""
    return [float(x) for x in s.split()]

# --------------------------------------------------------------------------------------
# Core parsing
# --------------------------------------------------------------------------------------

def parse_mphtxt(path: Path):
    """
    Parse COMSOL *.mphtxt and return (nodes, faces, tets, tet_tag, version).
    Relies on marker lines such as 'Mesh vertex coordinates', '# number of elements', and
    '# Geometric entity indice' to locate sections robustly.
    """
    lines = path.read_text().splitlines()

    num = []
    Pmesh = -1
    Pelem = []
    Pgeo  = []

    for i, line in enumerate(lines, start=1):
        if '#' in line and not line.startswith('#'):
            val = first_int_in_line(line)
            if val is not None:
                num.append(val)
        if 'Mesh vertex coordinates' in line or 'Mesh point coordinates' in line:
            Pmesh = i
        if '# number of elements' in line:
            Pelem.append(i)
        if '# Geometric entity indice' in line:
            Pgeo.append(i)

    if len(num) < 24 or Pmesh < 0 or len(Pelem) < 4 or len(Pgeo) < 4:
        raise RuntimeError("Unexpected mphtxt format; required markers not found.")

    version = num[3]
    nNodes  = num[5]
    nFaces  = num[18]
    nTets   = num[23]

    nodes = []
    for li in range(Pmesh + 1, Pmesh + nNodes + 1):
        vals = floats_from_line(lines[li-1])
        if len(vals) != 3:
            raise RuntimeError(f"Invalid node line at {li}")
        nodes.append(vals)

    faces = []
    for li in range(Pelem[2] + 2, Pelem[2] + nFaces + 2):
        idx = ints_from_line(lines[li-1])
        if len(idx) != 3:
            raise RuntimeError(f"Invalid face line at {li}")
        faces.append(tuple(idx))

    tets = []
    for li in range(Pelem[3] + 2, Pelem[3] + nTets + 2):
        idx = ints_from_line(lines[li-1])
        if len(idx) != 4:
            raise RuntimeError(f"Invalid tet line at {li}")
        tets.append(tuple(idx))

    tet_tag = []
    for li in range(Pgeo[3] + 1, Pgeo[3] + nTets + 1):
        idx = ints_from_line(lines[li-1])
        if not idx:
            raise RuntimeError(f"Invalid tet tag line at {li}")
        tet_tag.append(idx[0])

    return nodes, faces, tets, tet_tag, version

# --------------------------------------------------------------------------------------
# Gmsh v4 ASCII parser (minimal, nodes + triangles + tets + physical groups)
# --------------------------------------------------------------------------------------

def parse_msh_v4_ascii(path: Path):
    """
    Parse Gmsh MSH 4.x (ASCII) per spec:
      - $Nodes: blocks with (entityDim, entityTag), node tags then coordinates
      - $Elements: blocks with (entityDim, entityTag, elementType)
      - $Entities: map volume entityTag -> physicalTag(s) used as domain IDs
    We collect:
      nodes: [(x,y,z)]
      tris:  [(n0,n1,n2)] with their entityTag
      tets:  [(n0,n1,n2,n3)] and a tet_tag[] (domain id) derived from the volume 
             entity's first physical tag.
    Requires tetrahedra (3D). If no tetrahedra are present, conversion aborts.
    Reference: Gmsh file formats, $Nodes/$Elements/$Entities, elementType 2=triangle, 4=tetra.
    """
    lines = path.read_text().splitlines()
    i = 0
    N = len(lines)
    def seek(tag):
        nonlocal i
        while i < N and lines[i].strip() != tag:
            i += 1
        return i < N

    # Entities: record physical groups per volume entity
    vol_entity_to_phys: Dict[int, int] = {}
    if seek("$Entities"):
        i += 1
        header = lines[i].split(); i += 1
        nP, nC, nS, nV = map(int, header)
        # skip points and curves
        i += nP  # each point line
        i += nC  # each curve line
        # surfaces: skip their lines (variable length) robustly: read count, then consume that many blocks
        # However, for our needs only volumes' physical tags are required. We'll jump to volumes by scanning.
        # Conservative approach: iterate nS and consume the variable-tail per spec:
        for _ in range(nS):
            toks = lines[i].split(); i += 1
            # surfaceTag + bbox(6) + numPhysicalTags + ... + numBoundingCurves + ...
            numPT = int(toks[7])
            # consume physical tags (can spill to next lines if long -> they do not; remain on same line)
            # Reconstruct remaining counts:
            remain = toks[8 + numPT:]
            if remain:
                numBC = int(remain[0])
                # each BC is 1 int (possibly signed) on the same line already; nothing more to skip here
            # else: rare/no bounding curves line; nothing more to do
        # volumes:
        for _ in range(nV):
            toks = lines[i].split(); i += 1
            # volumeTag, bbox(6), numPhysicalTags, [phys...] , numBoundingSurfaces, [surf ...]
            volTag = int(toks[0])
            numPT = int(toks[7])
            phys = [int(x) for x in toks[8:8+numPT]] if numPT > 0 else []
            vol_entity_to_phys[volTag] = phys[0] if phys else 0
            # bounding surfaces count is toks[8+numPT], we don't need to advance further
        # $EndEntities
        while i < N and lines[i].strip() != "$EndEntities":
            i += 1
        if i < N: i += 1
    else:
        i = 0  # no Entities block (optional in 4.1); leave mapping empty -> defaults to 0

    # Nodes
    if not seek("$Nodes"):
        raise RuntimeError("Invalid .msh: missing $Nodes")
    i += 1
    hb = [int(x) for x in lines[i].split()]; i += 1
    numBlocks, numNodes, minTag, maxTag = hb
    tag_to_index: Dict[int, int] = {}
    nodes: List[Tuple[float,float,float]] = [None] * numNodes  # type: ignore
    # blocks
    tags_in_order: List[int] = []
    for _ in range(numBlocks):
        b = [int(x) for x in lines[i].split()]; i += 1
        entityDim, entityTag, parametric, numInBlock = b
        # node tags (could span multiple lines; in practice they are contiguous numbers)
        blk_tags: List[int] = []
        for _ in range(numInBlock):
            t = int(lines[i].strip()); i += 1
            blk_tags.append(t); tags_in_order.append(t)
        # coordinates for this block
        for k, t in enumerate(blk_tags):
            xyz = [float(x) for x in lines[i].split()[:3]]; i += 1
            idx = len(tag_to_index)
            tag_to_index[t] = idx
            nodes[idx] = (xyz[0], xyz[1], xyz[2])
            if parametric:
                # skip u,[v],[w] if present
                pass
    if lines[i].strip() != "$EndNodes":
        # advance to end if necessary
        while i < N and lines[i].strip() != "$EndNodes":
            i += 1
    i += 1

    # Elements
    if not seek("$Elements"):
        raise RuntimeError("Invalid .msh: missing $Elements")
    i += 1
    hb = [int(x) for x in lines[i].split()]; i += 1
    numEBlocks, numElems, emin, emax = hb
    tris: List[Tuple[int,int,int]] = []
    tri_ent: List[int] = []
    tets: List[Tuple[int,int,int,int]] = []
    tet_entity: List[int] = []
    ET_TRI, ET_TET = 2, 4
    for _ in range(numEBlocks):
        b = [int(x) for x in lines[i].split()]; i += 1
        entityDim, entityTag, etype, numInBlock = b
        for _ in range(numInBlock):
            parts = [int(x) for x in lines[i].split()]
            i += 1
            elTag = parts[0]
            nn = parts[1:]
            if etype == ET_TRI and len(nn) >= 3:
                a,bn,c = [tag_to_index[nn[0]], tag_to_index[nn[1]], tag_to_index[nn[2]]]
                tris.append((a,bn,c)); tri_ent.append(entityTag)
            elif etype == ET_TET and len(nn) >= 4:
                a,bn,c,d = [tag_to_index[nn[0]], tag_to_index[nn[1]], tag_to_index[nn[2]], tag_to_index[nn[3]]]
                tets.append((a,bn,c,d)); tet_entity.append(entityTag)
            else:
                # ignore other element types
                pass
    # drain to $EndElements
    while i < N and lines[i].strip() != "$EndElements":
        i += 1
    if i < N: i += 1

    # derive tet_tag (domain id) from volume entity -> first physical tag
    tet_tag: List[int] = []
    if tets:
        for e in tet_entity:
            tet_tag.append(vol_entity_to_phys.get(e, 0))
    else:
        raise RuntimeError("Gmsh .msh must be a 3D tetrahedral mesh (no 2D-only support).")            
    version = 4
    return nodes, tris, tets, tet_tag, tri_ent, version

# --------------------------------------------------------------------------------------
# Domain assignment
# --------------------------------------------------------------------------------------

def build_face_to_tets(tets: List[Tuple[int,int,int,int]]) -> Dict[Tuple[int,int,int], List[int]]:
    """Map each unique triangle (as sorted node triplet) to the list of adjacent tetra indices."""
    face2tets: Dict[Tuple[int,int,int], List[int]] = {}
    for ti, (a, b, c, d) in enumerate(tets):
        for f in ((a, b, c), (a, b, d), (a, c, d), (b, c, d)):
            key = tuple(sorted(f))
            face2tets.setdefault(key, []).append(ti)
    return face2tets

def assign_domains(nodes, faces, tets, tet_tag, periodic: bool, PP: List[float]):
    """
    Derive (front, back) domain IDs for each surface triangle.
    - Skip faces outside the periodic window when periodic=True.
    - Determine orientation by comparing face normal with vector to a neighboring tet center.
    """
    face2tets = build_face_to_tets(tets)
    domains = []
    for fi, (i, j, k) in enumerate(faces):
        if periodic and not test_edge_periodic(nodes, (i, j, k), PP):
            continue
        key = tuple(sorted((i, j, k)))
        adj = face2tets.get(key, [])
        ftag, btag = 0, 0
        if len(adj) >= 2:
            t0, t1 = adj[0], adj[1]
            ftag, btag = tet_tag[t0], tet_tag[t1]
            A, B, C = nodes[i], nodes[j], nodes[k]
            FC = tri_center(A, B, C)
            fn = face_normal(A, B, C)
            A0, B0, C0, D0 = nodes[tets[t0][0]], nodes[tets[t0][1]], nodes[tets[t0][2]], nodes[tets[t0][3]]
            TC = tet_center(A0, B0, C0, D0)
            if dot([TC[0]-FC[0], TC[1]-FC[1], TC[2]-FC[2]], fn) < 0:
                ftag, btag = btag, ftag
        elif len(adj) == 1:
            btag = tet_tag[adj[0]]
        domains.append((fi, ftag, btag))
    return domains

# --------------------------------------------------------------------------------------
# Writer
# --------------------------------------------------------------------------------------

def write_mesh_mesh(out_path: Path, nodes, faces, domains):
    """
    Write HELIOS mesh format:
    Nnodes
    x y z (per node)
    Nfaces
    i j k  front back   (1-based indices; enforce front<=back by swapping i<->j if needed)
    """
    with out_path.open('w') as f:
        f.write(f"{len(nodes)}\n")
        for x, y, z in nodes:
            f.write(f"{x:.8f}\t{y:.8f}\t{z:.8f}\n")
        f.write(f"{len(domains)}\n")
        for face_idx, front, back in domains:
            i, j, k = faces[face_idx]
            if front <= back:
                f.write(f"{i+1}\t{j+1}\t{k+1}\t{front}\t{back}\n")
            else:
                f.write(f"{j+1}\t{i+1}\t{k+1}\t{back}\t{front}\n")

# --------------------------------------------------------------------------------------
# New helpers: project-aware discovery
# --------------------------------------------------------------------------------------

def clean_lines(path: Path) -> List[str]:
    """Read config.txt, removing blank lines and full-line comments."""
    out = []
    with open(path, 'r') as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith('#'):
                continue
            out.append(s)
    return out

def detect_periodic_from_config(cfg: Path):
    """
    Heuristically detect periodic parameters from config.txt.
    If a line with exactly 5 float-like tokens exists -> treat as (px,py,cx,cy,zcut).
    """
    try:
        lines = clean_lines(cfg)
    except Exception:
        return 0, []
    for s in lines:
        toks = s.replace(',', ' ').split()
        if len(toks) == 5:
            try:
                vals = list(map(float, toks))
                return 1, vals
            except ValueError:
                continue
    return 0, []

def find_mesh_file(input_dir: Path) -> Path:
    """
    Find *.mphtxt OR *.msh in input_dir.
    Prioritizes 'mesh.mphtxt' or 'mesh.msh' if multiple exist.
    """
    cands = sorted(list(input_dir.glob('*.mphtxt')) + list(input_dir.glob('*.msh')))
    if not cands:
        raise FileNotFoundError(f"No .mphtxt or .msh found in {input_dir}")
    
    # Prefer a file already named mesh.*
    for c in cands:
        if c.stem.lower() == 'mesh':
                return c
    return cands[0]    

# --------------------------------------------------------------------------------------
# Conversion wrapper
# --------------------------------------------------------------------------------------

def convert_one(src: Path, periodic_flag: int, PP: List[float], out_path: Path = None):
    """
    End-to-end conversion: parse .mphtxt or .msh -> assign domains 
    (with optional periodic window) -> write mesh.mesh.
    """
    nodes = []; faces = []; tets = []; tet_tag = []; version = 0
    tri_entity: Optional[List[int]] = None
    if src.suffix.lower() == ".mphtxt":
        nodes, faces, tets, tet_tag, version = parse_mphtxt(src)
        if version != 4:
            print("Warning: expected COMSOL >= 5.x mphtxt; proceeding anyway.")
        print(f"Data acquired with: {len(nodes)} nodes, {len(faces)} faces and {len(tets)} tetras")
        domains = assign_domains(nodes, faces, tets, tet_tag, periodic=bool(periodic_flag), PP=PP)
    else:
        nodes, faces, tets, tet_tag, tri_entity, version = parse_msh_v4_ascii(src)
        print(f"Gmsh mesh detected: nodes={len(nodes)} tris={len(faces)} tets={len(tets)}")
        domains = assign_domains(nodes, faces, tets, tet_tag, periodic=bool(periodic_flag), PP=PP)
    print("Assigning domains: 100 %")
    print("Triangle data sorted.")
    if out_path is None:
        out_path = Path('mesh.mesh')
    out_path.parent.mkdir(parents=True, exist_ok=True)
    print("Writing converted mesh...")
    write_mesh_mesh(out_path, nodes, faces, domains)
    max_tag = max(tet_tag) if tet_tag else 0
    print(f"Mesh successfully converted with {len(domains)} triangles and {max_tag+1} domains!")
    print(f"Output: {out_path.resolve()}")

# --------------------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------------------

def main():
    args = sys.argv[1:]

    # ------------------------------------------------------------------
    # New short form: meshconvert.py SimName
    # ------------------------------------------------------------------
    if len(args) == 1 and not args[0].endswith('.mphtxt') and args[0] not in ('-o', '--output'):
        sim_name = args[0]
        in_dir  = (SIM_DATA_ROOT / sim_name).resolve()
        out_dir = (SIM_RES_ROOT  / sim_name).resolve()

        if args[0] in ('-h', '--help'):
            print('Usage (new):   meshconvert.py SimName')
            print('Usage (legacy):')
            print('  Isolated: meshconvert.py [-o OUTPATH] [MeshFilePath.mphtxt]')
            print('  Periodic: meshconvert.py [-o OUTPATH] [MeshFilePath.mphtxt] px py cx cy zcut')
            sys.exit(0)

        cfg = in_dir / 'config.txt'
        if not cfg.is_file():
            print(f"Error: config.txt not found in {in_dir}")
            sys.exit(2)

        # Detect (periodic_flag, PP) from config
        periodic_flag, PP = detect_periodic_from_config(cfg)

        # Discover mesh, copy to sim_res/SimName/mesh.<ext>
        src_mesh = find_mesh_file(in_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        dst_mesh = out_dir / src_mesh.name
        shutil.copy2(src_mesh, dst_mesh)
        print(f"Copied mesh: {src_mesh.name} -> {dst_mesh}")

        # Convert to sim_res/SimName/mesh.mesh
        out_mesh = out_dir / 'mesh.mesh'
        print("Mesh found.")
        convert_one(dst_mesh, periodic_flag, PP, out_mesh)
        return

    # ------------------------------------------------------------------
    # Legacy forms (unchanged)
    # ------------------------------------------------------------------
    # Optional -o/--output
    out_path = None
    if "-o" in args:
        i = args.index("-o")
        if i + 1 >= len(args):
            print("Error: -o requires a path")
            sys.exit(2)
        out_path = Path(args[i + 1])
        del args[i:i + 2]
    elif "--output" in args:
        i = args.index("--output")
        if i + 1 >= len(args):
            print("Error: --output requires a path")
            sys.exit(2)
        out_path = Path(args[i + 1])
        del args[i:i + 2]

    # Periodic: 5 numeric args (no file)
    if len(args) == 5:
        try:
            px, py, cx, cy, zcut = map(float, args)
        except ValueError:
            print("Error: periodic parameters must be numeric: px py cx cy zcut")
            sys.exit(2)
        PP = [px, py, cx, cy, zcut]
        periodic_flag = 1
        mphtxts = [f for f in os.listdir('.') if os.path.isfile(f) and f.endswith('.mphtxt')]
        if not mphtxts:
            print("No .mphtxt found in current directory.")
            sys.exit(1)
        src = Path(mphtxts[0])
        print("Mesh found.")
        convert_one(src, periodic_flag, PP, out_path)
        return

    # Periodic: explicit file + 5 numbers
    if len(args) == 6:
        src = Path(args[0])
        if not src.is_file():
            print("File does not exist:", src.resolve()); sys.exit(1)
        try:
            px, py, cx, cy, zcut = map(float, args[1:])
        except ValueError:
            print("Error: periodic parameters must be numeric: px py cx cy zcut"); sys.exit(2)
        PP = [px, py, cx, cy, zcut]
        print("Mesh found.")
        convert_one(src, 1, PP, out_path)
        return

    # Isolated/Direct file: 0 or 1 arg (file optional). Supports .mphtxt or .msh
    if len(args) < 2:
        periodic_flag = 0
        PP = []
        if len(args) == 1:
            src = Path(args[0])
            if not src.is_file():
                print("File does not exist:", src.resolve()); sys.exit(1)
        else:
            # try common names
            cands = [f for f in os.listdir('.') if os.path.isfile(f) and (f.endswith('.mphtxt') or f.endswith('.msh'))]
            if not cands:
                print("No .mphtxt/.msh found in current directory."); sys.exit(1)
            src = Path(cands[0])
        print("Mesh found.")
        convert_one(src, periodic_flag, PP, out_path)
        return

if __name__ == "__main__":
    main()
