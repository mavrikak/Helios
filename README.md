![alt text](third_party/logo_README.png)
# HELIOS - HomogEneous and Layered medIa Optical Scattering

HELIOS is a C++ / Python toolkit for simulating electromagnetic scattering from objects embedded either in **homogeneous** media (isolated objects and periodic structures) or **layered** media, using the **surface integral equation (SIE)** method. It provides a modern, scriptable workflow around two compiled applications (`SIENano` and `SIENanoPP`) and a set of Python utilities for job generation, mesh conversion, post-processing, and visualization.

---

## 🗂️ Table of contents

- [Citations](#-citations)
- [Highlights](#-highlights)
- [Repository layout](#-repository-layout)
- [Dependencies](#-dependencies)
- [Quick start](#-quick-start)
- [Input data and configuration](#-input-data-and-configuration)
- [Running simulations](#-running-simulations)
- [Visualization](#-visualization)
- [Graphical user interface](#-graphical-user-interface)
- [Logging, progress, and exit codes](#-logging-progress-and-exit-codes)

---

## 📝 Citations

- P.S. Mavrikakis, O.J.F. Martin, arXiv (2026). https://doi.org/10.48550/arXiv.2602.23097.
- A.M. Kern, O.J.F. Martin, J. Opt. Soc. Am. A 26 (2009) 732–740. https://doi.org/10.1364/JOSAA.26.000732.
- B. Gallinet, A.M. Kern, O.J.F. Martin, J. Opt. Soc. Am. A Opt. Image Sci. Vis. 27 (2010) 2261–2271. https://doi.org/10.1364/JOSAA.27.002261.

---

## ✨ Highlights

- **Three modes**: 
  - isolated scatterer in homogeneous background, 
  - periodic structure in homogeneous background, and 
  - isolated scatterer in layered background.
- **Scriptable orchestration** with `run_sie.py`: prepare jobs, run solvers, and post-process results from the command line.
- **Robust job generation** from a compact `config.txt` using `pytools/jobwriter.py`, including layered stacks and tabulated materials.
- **COMSOL mesh imports**: convert `.mphtxt` to the compact `.mesh` format with `pytools/meshconvert.py`, including optional periodic cropping.
- **Near-field & cross-section visualization** with `pytools/visualization.py` in 2D, 3D, and 1D line-cut modes.
- Designed to work with **C++11 builds** of `apps/SIENano` and `apps/SIENanoPP`.

---

## 🧭 Repository layout

```bash
Helios/
├── docs/                  # C++ code documentation folder
├── include/               # header files and blitz++ files
├── materials/             # ε(λ) tables for named materials
├── pytools/               # python orchestration & utilities
│   ├── drude_model.py         # Drude model for tabulated materials
│   ├── jobwriter.py           # turns sim_data/<sim>/config.txt into jobs/job.N.xml
│   ├── meshconvert.py         # COMSOL .mphtxt → .mesh
│   ├── plot_mesh.py           # plot simulation mesh
│   └── visualization.py       # plot cross sections and near-fields
├── sim_data/              # your data inputs (per simulation)
│   └── <SimName>/           # simulation data folder(s)
│       ├── config.txt         # compact configuration for jobwriter.py
│       └── mesh.mphtxt        # COMSOL text mesh (.mphtxt)
├── sim_res/               # data outputs (auto-created per simulation)
│   └── <SimName>/           # simulation results folder(s)
│       ├── jobs/              # folder where job.N.xml files are found
│       ├── logs/              # one logfiles per run (different for solver/post-processor)
│       ├── out/               # solver & post-processor data outputs
│       │  ├─ csc/             # cross sections' data
│       │  ├─ fields/          # .ein/.esc/.hin/.hsc fields' data organized by points set
│       │  └─ media/           # location for PNG plots
│       ├── points/            # generated sampling points sets (.pos)
│       ├── lambdalist.txt     # λ[nm] mapping for plots
│       └── mesh.mesh          # converted surface mesh
├── src/                   # C++ source code
├── third_party/           # third-party fortran files
│   ├── amos/                   # amos library files
│   ├── legacy_f                # congugate gradient solver
│   ├── quadpack                # quadpack routines
│   ├── regridpack              # regridpack routines
│   └── erfc.bin                # erfc lookup table for periodic jobs
├── dist_run_helios.py     # parallel job launcher for distributed runs
├── environment.yml        # conda env (see below)
├── helios_gui.py          # Python GUI
├── makefile               # compilation instructions
├── README.md              # this file
└── run_sie.py             # orchestration: prepare → solve → generate points → post-process
```

---

## ⚙️ Dependencies

### System Overview

The project consists of:
- **C++11 core code** (compiled via `makefile`),
- **Fortran numerical kernels** (AMOS, QUADPACK, REGRIDPACK),
- **Python tools** (for job management, post-processing, and visualization),
- **Optional bash launchers** (for distributed runs).

>All builds have been tested on Linux using `conda` compilers and OpenBLAS.

### 🧰 Core Build Dependencies

| Category | Package | Purpose |
|-----------|----------|----------|
| **Compiler (C++)** | `gxx_linux-64` | C++11 compiler used for all main code. |
| **Compiler (Fortran)** | `gfortran_linux-64` | Fortran compiler for numerical submodules. |
| **Math Library** | `openblas` | Provides BLAS and LAPACK routines. |
| **C++ Array Library** | `blitz++` | High-performance multi-dimensional arrays used in the solver. |
| **Threads** | `pthread` *(system)* | Required by OpenBLAS. |

**Fortran sources included in the project:**
- `third_party/amos/*.f` — Bessel and Hankel functions,
- `third_party/quadpack/*.f` — Gauss-Kronrod integration,
- `third_party/regridpack/*.f90` — interpolation routines,
- `third_party/legacy_f/cg.f` — conjugate-gradient solver.

>These are compiled automatically and do **not** require separate installation.

### 🐍 Python Dependencies

Used by the helper tools under `pytools/` and for running the main `run_sie.py` script and `helios_gui.py`.

| Package | Purpose |
|----------|----------|
| `python` ≥ 3.10 | Base interpreter. |
| `numpy` | Numerical processing and mesh handling. |
| `scipy` | Numerical routines. |
| `matplotlib` | Plotting and visualization. |
| `PySide6` | Qt bindings. |
| `tqdm` | Progress bars. |
| `miepython` | Mie theory. |

### 📖 Documentation Dependencies

The documentation is already provided, but the following packages are required to regenerate it with **Doxygen**.

| Package | Purpose |
|----------|----------|
| `doxygen` | Generates HTML and PDF documentation from C++ comments. |
| `graphviz` | Renders dependency and call graphs (used by Doxygen). |

### 🧱 Recommended Conda Environment

The project uses [Miniconda](https://www.anaconda.com/docs/getting-started/miniconda/install#linux-2) to manage dependencies. The file `environment.yml` is provided to produce the proper conda environment. After installing [Miniconda](https://www.anaconda.com/docs/getting-started/miniconda/install#linux-2) for linux, you can create and activate the environment with the following commands:

```bash
cd Helios/
conda env create -f environment.yml
conda activate Helios
```

---

## 🚀 Quick start

***Linux users can clone the repository and install HELIOS with the following steps. For Windows users, the [Windows Subsystem for Linux (WSL)](https://learn.microsoft.com/en-us/windows/wsl/install) is recommended.***

1. Install git and clone the repository:

```bash
sudo apt install git git-lfs && git lfs install
git clone https://github.com/mavrikak/Helios.git
```

2. Create a conda environment with [Miniconda](https://www.anaconda.com/docs/getting-started/miniconda/install#linux-2) and activate it:

```bash
cd Helios
conda env create -f environment.yml
conda activate Helios
```

3. Build the project:

```bash
mkdir apps
mkdir build
make clean && make -j
```

4. Place inputs for a simulation under `sim_data/<SimName>/`:

    - `config.txt` — simulation configuration (check examples in `sim_data/` and next section for details)
    - `mesh.mphtxt` — COMSOL mesh text file

5. [Optional] Rebuild the documentation:

```bash
cd docs
bash build.sh
```

6. [Optional] Generate `erfc(z)` lookup table:

```bash
make table
python run_sie.py make-table erfc.bin --maxRe 10.0 --maxIm 10.0 --incRe 0.1 --incIm 0.1
```

---

## 🧩 Input data and configuration

Each simulation lives under `sim_data/<SimName>/` and provides a `config.txt`. From this compact text, `pytools/jobwriter.py` produces the solver job XMLs. The configuration supports:

- Wavelength sweep: Specify one value or a `start stop step` value range.

- Domain description: tokens for each domain; plain materials can be either literals (real/complex numbers) or names of tabulated materials found in `materials/`.

- Incident domain: the index of the domain where the illuminating field lives.

- Polarization: `s, p, TE, TM, LCP, RCP`.

- Incidence angles: $\theta$ and $\phi$ in degrees; specify one value or a `start stop step` value range.

### Layered media blocks

Add two extra lines to list, per layered stack:

- Materials per stack (comma-separated stacks, space-separated within a stack): `1 SiO2 1, Ag Au`

- Interface z-positions per stack (same structure): `0 150, 20`

- jobwriter.py expands each `LayeredM#` token in your domain line into a `<domainLayers>` section with `<layer>` entries (properties at the current wavelength) and `<interface>` tags.

### Materials

Material tokens in `config.txt` can be:

- Literals: $\varepsilon_r = \alpha + i \beta$ (`2.25, 2.25+0.1i, 2.25+0.1j, (2.25, 0.1), [2.25,0.1], {2.25,0.1}`).

- Tabulated files in `materials/<name>.txt` with the format: $ \quad \lambda \, [nm] \quad \Re(\varepsilon_r) \quad \Im(\varepsilon_r) \quad$.

>Linear interpolation fills intermediate wavelengths; endpoints are linearly extrapolated.

### Periodic structure blocks

Add one extra line to list:

- `(px, py)`: Unit cell dimensions in the x and y directions, respectively.
- `(cx, cy)`: Unit cell centers in the x and y directions, respectively.
- `z`: Cut-off z position for the lower semi-infinite part of the unit cell.

---

## 🖱 Graphical user interface

The `helios_gui.py` script can be used to open the GUI and interactively prepare simulations, run them (both locally and distributed), and post-process the results. Additionally, it can be used to visualize the mesh and results of a simulation. All the user has to do is run `python3 helios_gui.py` in a terminal.

---

## 📜 Logging, progress, and exit codes

- Every external run is tee-logged under `sim_res/<sim>/logs/` with concise TTY progress when interactive.

- Non-zero exit codes indicate user-facing errors (bad paths/flags/inputs).

- On non-interactive terminals, spinners are suppressed to keep logs clean.

---