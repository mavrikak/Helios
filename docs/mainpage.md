# HELIOS - HomogEneous and Layered medIa Optical Scattering

# Overview

**HELIOS (HomogEneous and Layered medIa Optical Scattering)** is a C++ codebase for accurate and efficient simulations of electromagnetic scattering from objects embedded in homogeneous or layered media. It implements the surface integral equation (SIE) method with the Poggio–Miller–Chang–Harrington–Wu-Tsai (PMCHWT) formulation and Rao–Wilton–Glisson (RWG) as both test and basis functions in a Galerkin scheme. Additionally, it supports periodic boundary conditions for homogeneous backgrounds, and provides robust machinery for handling singular and near-singular integrals that arise. The software emphasizes clarity, numerical stability, and extensibility: geometry and meshing are cleanly separated from physics (Green’s functions, material models, incident fields) and from solvers and workflows.

At a high level, the program takes a surface mesh describing one or more scatterers, assigns domain indices to each triangle (front/back), and assembles the linear system corresponding to the PMCHWT SIE formulation for the unknown equivalent surface currents. From these currents, it evaluates near/far fields, scattering parameters, and other derived quantities, in either a single homogeneous background or a stratified stack (multi-layer), where Sommerfeld integrals are required.

See also: \ref modules

# Key capabilities

- **Formulations**

    - PMCHWT for penetrable scatterers.

    - Clean base class (SIEForm) that factors common assembly patterns and allows new formulations to be plugged in.

- **Basis and geometry**

    - RWG basis on triangular surface meshes.

    - Geometry primitives (Triangle, SurfaceMesh) cache normals, areas, centroids, self-integrals, and adjacency information.

- **Green’s functions**

    - Homogeneous 3D Green’s functions for free space and homogeneous media.

    - Periodic homogeneous 3D Green’s functions for 2D periodicities.

    - Layered-media Green’s functions with Sommerfeld integral evaluation and stabilization; utilities for interface indexing, layer properties, and wavenumber bookkeeping.

- **Incident fields**

    - Plane wave and dipole sources (with extensible design for additional sources).

    - Layer-aware evaluation routines for E/H in stratified stacks.

- **Quadrature and singularity treatment**

    - Standard Gauss rules and polar quadrature (TriQuadPol) for near-singular panels.

    - Analytic singularity treatment pieces (SingSub, SingCancellation).

- **Layered media workflow**

    - Utilities to translate physical layer stacks into the frequency-domain quantities used by Green’s functions and incident fields.

    - Efficient reuse of precomputed quantities across many evaluation points.

- **I/O and job orchestration**

    - Mesh I/O for custom .mesh, COMSOL .mphtxt, and Gmsh .msh formats.
    
    - Simple "job parser" to script conversions (e.g., plane cuts, domain relabeling) and to keep assembly/runs reproducible.

# Architecture in brief

The code is divided into loosely coupled modules:

- **Mesh & Geometry:** Management of nodes, triangles, edges, domain indices, normals/areas, self-integrals, and neighbor queries. This module owns nothing about physics—just geometry and indexing.

- **Green’s Functions:** Implemention of scalar/tensor Green’s functions and related transforms. They are the "physics kernels" used during assembly and field evaluation.

- **Formulations:** Orchestrate element–by–element contributions, call the appropriate Green’s functions, and assemble the global linear system (left-hand side A matrix and right-hand side vector b).

- **Incident Fields:** Provides E/H evaluations at arbitrary points and their integration with RWG functions for the excitation vector.

- **Quadrature & Singularity Handling:** Gauss rules, polar decomposition around a query point, and analytic subtractions, cancellations and kernel refinements.

- **Driver: `SimJob` ties everything together:** parsing inputs, constructing the mesh, selecting the formulation, setting up sources, building/solving the linear system, and writing outputs.

This separation keeps modules testable in isolation and makes it straightforward to add new formulations, sources, or Green’s-function back-ends without disturbing the rest of the code significantly.

# Numerical notes

- **Accuracy near singularities:** Close-panel interactions are handled using cancellation methods and analytic subtractions of the kernel’s singular part and by specialized polar quadratures. This preserves accuracy as observation points approach element edges or vertices.

- **Layered integrals:** Sommerfeld integrals are treated with stable contour choices and asymptotic safeguards. Layer parameters (wavenumbers, impedances) are cached and reused across elements to reduce overhead in multi-query scenarios.

- **Orientation and domains:** Each triangle carries explicit front/back domain indices; reference normals point to the front domain. This convention is enforced consistently across basis assembly and singular terms.

# Typical workflow

**Load mesh:** Read a .mesh, .mphtxt, or .msh file into SurfaceMesh; verify domain indices and orientations; optionally run mesh conversions (plane trimming, index remapping).

**Choose physics:** Select homogeneous (isolated structure or periodic simulation) or layered medium; set material parameters for each domain or layer; configure Green’s-function back-end.

**Set incident field:** Create a PlaneWave or Dipole (or Gaussian for homogeneous backgrounds) with direction, polarization, and frequency; for layered cases, provide the layer stack description.

**Assemble system**
Pick an SIE formulation (currently only PMCHWT); assemble the system matrix and excitation vector using the mesh, basis, Green’s functions, and incident field; the assembly will invoke singular/near-singular handling as needed.

**Solve**
Use your your preffered direct or iterative solver (external to this core) to obtain surface current coefficients.

**Post-process**
Evaluate near/far fields, cross sections, or custom functionals by re-using the Green’s functions and basis projections.

# Inputs and outputs

**Inputs:** mesh files, layer stacks (permittivity/permeability vs. depth), incident field parameters, and solver options.

**Outputs:** current coefficients, E/H near-/far-fields, and optional CSV dumps for visualization or downstream analysis.

# Extensibility

**New incident fields:** take advantage of the IncidentField interface and add new illuminations as new children classes.

**New formulations:** derive from SIEForm, reusing quadrature and Green’s-function hooks.

**Alternate Green’s functions:** take advantage of the GreenF interface and add new Green's function implementations as new children classes.

**Additional quadratures:** plug into the same evaluation sites by following the existing patterns.

# Applications

**Executable files:** All the application files of the HELIOS project are located in folder 'apps/'. These are the executable files for the calculation of current coefficients and cross sections (SIENano), near-/far-fields at a set of points (SIENanoPP), and Casimir forces (Casimir).

# Documentation

- The codebase is documented with Doxygen comments at file, class, and function levels.

- The documentation site (generated into docs/build/html/) includes browsable source, call/include graphs (if Graphviz is enabled), module pages, and a class index.

- Group tags (`\ingroup`) organize classes into logical modules so readers can navigate from a conceptual overview to implementation details.

# Limitations and future work

- The current focus is frequency-domain scattering with RWG bases on triangular meshes; volumetric effects and transient behavior are out of scope.

- Solver infrastructure (preconditioners, fast methods) is deliberately left open so you can integrate your favorite linear algebra stack!
