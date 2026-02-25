\page groups Groups

This page links to the main modules that organize the codebase.  
Use the left navigation (**Topics → Modules**) to browse each module, or click below:

- \ref mesh   "Mesh & Geometry"
- \ref domain "Domains"
- \ref greenFunction "Green's Functions"
- \ref formulation   "SIE Formulations"
- \ref incidentField    "Incident Fields"
- \ref quadrature   "Quadrature & Singularity Handling"
- \ref driver "Driver & Jobs"

---

\defgroup modules Modules
High-level modules that organize the code base.

\defgroup mesh Mesh & Geometry
\ingroup modules
@{
Surface mesh and geometric primitives.
@}

\defgroup domain Simulation domains
\ingroup modules
@{
Homogeneous and layered domains.
@}

\defgroup greenFunction Green's Functions
\ingroup modules
@{
Homogeneous and layered Green’s back-ends and Sommerfeld integration.
@}

\defgroup formulation SIE Formulations
\ingroup modules
@{
System assembly and gradients.
@}

\defgroup incidentField Incident Fields
\ingroup modules
@{
Incident sources and their integrations.
@}

\defgroup quadrature Quadrature & Singularity Handling
\ingroup modules
@{
Special quadratures and analytic subtractions.
@}

\defgroup driver Driver & Jobs
\ingroup modules
@{
Orchestration, file I/O, and batch runs.
@}
