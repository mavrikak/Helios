/**
 * @file SimJob.h
 * @brief Orchestrates a full simulation run: parse job/subjobs, build the
 *        linear system, solve it, and optionally post-process fields.
 *
 * One SimJob corresponds to one mesh and wavelength setup, with possibly
 * multiple subjobs (different incident fields) reusing the same system matrix.
 */

// Only include this once during compiling
#ifndef SIMJOB_H
#define SIMJOB_H

#define EFIELD 1   ///< Output electric field E
#define HFIELD 2   ///< Output magnetic field H
#define DFIELD 4   ///< Output electric displacement @f$\mathbf{D} = \varepsilon \mathbf{E}@f$
#define BFIELD 8   ///< Output magnetic flux density @f$\mathbf{B} = \mu \mathbf{H}@f$
#define POYVEC 16  ///< Output Poynting vectors (inc/sca/ext)
#define FORCE 32   ///< (reserved) force output
#define POWER 64   ///< Output total power crossing domain boundaries

#include <blitz/array.h>
#include <mutex>
#include <string>
#include <vector>
#include "Domain.h"
#include "JobParser.h"
#include "SIEForm.h"

/// \ingroup driver
/**
 * @class SimJob
 * @brief High-level driver for SIE simulations and field evaluation.
 *
 * Workflow:
 *  1) LoadFromFile() parses the job file: mesh, materials, sources, etc.
 *  2) CreateLSEMatrix() allocates the A matrix; domains fill it via the
 *     formulations (PMCHWT).
 *  3) For each subjob: CreateLSEVec() and CreateSolVec(), fill excitation
 *     vectors, then SolveLSE() to get the solution current vector(s).
 *  4) Optionally FieldEval() to sample fields or power/cross sections.
 */
class SimJob {
 protected:
  // --- Model state shared across a job ---
  std::vector<Domain*> domains;                     ///< All problem domains
  SurfaceMesh* mesh;                                ///< Shared surface mesh
  blitz::Array<dcmplx, 2>* lseMatrix;               ///< System matrix A (2Ne x 2Ne)
  blitz::Array<dcmplx, 2>* lseGradientMatrix[3];    ///< @f$\partial Z / \partial p_i@f$ matrices (Casimir)
  std::vector<blitz::Array<dcmplx, 1>*> lseVectors; ///< RHS vectors b
  std::vector<blitz::Array<dcmplx, 1>*> solVectors; ///< Solution vectors x
  std::vector<std::string> labels;                  ///< Labels per (sub)job
  dcmplx vacuumWavelength;                          ///< @f$\lambda_0@f$ used to compute @f$\omega@f$
  rvec blochVector;                                 ///< Bloch vector (periodic jobs)
  std::vector<rvec> latticeVectors;                 ///< Lattice vectors (periodic)
  JobParser jobParser;                              ///< XML-like job reader
  std::vector<IncidentField*> subjobIncidentFields; ///< Temp storage during load
  bool isCasimirJob;                                ///< Determine if this is a Casimir job
  static int threads;                               ///< Worker threads for field eval

  // ---- Parsers used inside LoadFromFile() ----
  /**
   * @brief Parse and load simulation label from @f$<label>...</label>@f$ section.
   * @return 0 on success; non-zero on parse error.
   */
  int LoadLabel();

  /**
   * @brief Parse and load mesh data from @f$<mesh>...</mesh>@f$ (non-periodic case).
   * @return 0 on success; non-zero on parse error.
   */
  int LoadMesh();

  /**
   * @brief Parse and load periodic mesh from @f$<mesh>...</mesh>@f$ (applies symmetry operations).
   * @return 0 on success; non-zero on parse error.
   */
  int LoadMeshPer();

  /**
   * @brief Parse and load wavelength from @f$<wavelength>...</wavelength>@f$.
   * @return 0 on success; non-zero on parse error.
   */
  int LoadWavelength();

  /**
   * @brief Parse and load homogeneous domain definition from @f$<domain>...</domain>@f$.
   * @return 0 on success; non-zero on parse error.
   */
  int LoadDomain();

  /**
   * @brief Parse and load layered-media domain from @f$<domainLayers>...</domainLayers>@f$.
   * @return 0 on success; non-zero on parse error.
   */
  int LoadDomainLayered();

  /**
   * @brief Parse and load periodic domain definition from @f$<periodicdomain>...</periodicdomain>@f$.
   * @return 0 on success; non-zero on parse error.
   */
  int LoadDomainPer();

  /**
   * @brief Parse and load plane-wave excitation from @f$<planewave>...</planewave>@f$.
   * @return 0 on success; non-zero on parse error.
   */
  int LoadPlanewave();

  /**
   * @brief Parse and load dipole excitation from @f$<dipole>...</dipole>@f$.
   * @return 0 on success; non-zero on parse error.
   */
  int LoadDipole();

  /**
   * @brief Parse and load Gaussian-beam excitation from @f$<gaussian>...</gaussian>@f$.
   * @return 0 on success; non-zero on parse error.
   */
  int LoadGaussian();

  /**
   * @brief Parse and load periodicity parameters from @f$<periodic>...</periodic>@f$.
   * @return 0 on success; non-zero on parse error.
   */
  int LoadPeriodicity();

  LookupTableBin* t;        ///< Look-up table for the periodic Green's function computation
  bool thereIsALookupTable; ///< True if there is a lookup table

 public:
  // ---------------- Lifecycle ----------------
  /** @brief Construct an empty job (no mesh/domains yet). */
  SimJob();
  /** @brief Destructor: frees owned domains and LUT if allocated. */
  virtual ~SimJob();
  
  // ---------------- Job I/O ----------------
  /**
   * @brief Parse a full @f$<job>...</job>@f$ block from the already-open JobParser.
   * @return 0 on success; non-zero on parse/consistency error.
   */
  int LoadFromFile();
  
  /**
   * @brief Parse one @f$<subjob>...</subjob>@f$ (incident specification + label).
   * @return 0 on success; non-zero if malformed or missing label.
   */
  int LoadSubjobFromFile();
  
  // ---------------- Linear system construction ----------------
  /**
   * @brief Allocate and zero the system matrix A (and gradients if enabled).
   *        Matrix size is 2Ne x 2Ne, where Ne is mesh edge count.
   * @return 0 on success.
   */
  int CreateLSEMatrix();
  
  /**
   * @brief Allocate and zero a new RHS vector b, and append to lseVectors.
   * @return 0 on success.
   */
  int CreateLSEVec();

  /**
   * @brief Allocate and zero a new solution vector x, and append to solVectors.
   * @return 0 on success.
   */
  int CreateSolVec();
  
  // ---------------- Solvers ----------------
  /**
   * @brief Solve the current linear system with the chosen solver.
   * @param solverType 0: direct; 1: cgsqr; 2: LU (default).
   * @return 0 on success; non-zero on failure.
   */
  int SolveLSE(int solverType);
  
  /**
   * @brief Write a 1D Blitz vector in text format.
   * @tparam T Element type (e.g., dcmplx).
   * @param fileName Output filename.
   * @param vector   Vector to serialize.
   * @return 0 on success; non-zero on I/O error.
   */
  template <class T>
  int WriteVectorToFile(std::string fileName, blitz::Array<T, 1>* vector);
  
  /**
   * @brief Load a previously saved solution vector into solVectors.back().
   * @param solFile Input filename (.sol).
   * @return 0 on success; non-zero on I/O/format error.
   */
  int LoadSolFromFile(std::string solFile);

  // ---------------- Orchestration ----------------
  /**
   * @brief High-level entry: run a simulation for a given job file.
   * @param jobFile    Path to job file.
   * @param solverType Solver selection (as in SolveLSE()).
   * @param outPrefix  Prefix for output filenames.
   * @return 0 on success; non-zero on any failure.
   */
  int Simulate(std::string jobFile, int solverType, std::string outPrefix);
  
  // ---------------- Point-wise queries ----------------
  /**
   * @brief Find which domain contains a point.
   * @param pos Cartesian position.
   * @return Pointer to the containing Domain (background if none), never null.
   */
  Domain* ContainingDomain(rvec pos);
  
  /**
   * @brief Incident and secondary @f$\mathbf{E}@f$ at pos (domain auto-detected).
   * @param pos Cartesian position.
   * @return @f$\{\mathbf{E}_\mathrm{inc}(pos), \mathbf{E}_\mathrm{scat}(pos)\}@f$.
   */
  std::pair<cvec, cvec> EField(rvec pos);

  /**
   * @brief Incident and secondary @f$\mathbf{E}@f$ at pos with known containing domain.
   * @param pos Cartesian position.
   * @param contDomain Containing domain.
   * @return @f$\{\mathbf{E}_\mathrm{inc}(pos), \mathbf{E}_\mathrm{scat}(pos)\}@f$.
   */
  std::pair<cvec, cvec> EField(rvec pos, Domain* contDomain);
  
  /**
   * @brief Incident and secondary @f$\mathbf{H}@f$ at pos (domain auto-detected).
   * @param pos Cartesian position.
   * @return @f$\{\mathbf{H}_\mathrm{inc}(pos), \mathbf{H}_\mathrm{scat}(pos)\}@f$.
   */
  std::pair<cvec, cvec> HField(rvec pos);

  /**
   * @brief Incident and secondary @f$\mathbf{H}@f$ at pos with known containing domain.
   * @param pos Cartesian position.
   * @param contDomain Containing domain.
   * @return @f$\{\mathbf{H}_\mathrm{inc}(pos), \mathbf{H}_\mathrm{scat}(pos)\}@f$.
   */
  std::pair<cvec, cvec> HField(rvec pos, Domain* contDomain);

  /**
   * @brief Return both @f$(\mathbf{E}@f$,@f$\mathbf{H})@f$ incident and secondary pairs at pos.
   * @param pos Cartesian position.
   * @param contDomain Containing domain.
   * @return @f$\{\{\mathbf{E}_\mathrm{inc}(pos),  \mathbf{H}_\mathrm{inc}(pos)\}, 
   *              \{\mathbf{E}_\mathrm{scat}(pos), \mathbf{H}_\mathrm{scat}(pos)\}\}@f$.
   */
  std::pair<std::pair<cvec, cvec>, std::pair<cvec, cvec> > EandHField(
      rvec pos, Domain* contDomain);
  
  // ---------------- Diagnostics ----------------
  /** 
   * @brief Net time-averaged power crossing the domain's surface. 
   * @param domain Domain over which to compute @f$ P @f$.
   * @return Net power.
   * @note Currently only implemented for homogeneous domains.
   */
  double Power(Domain* domain);
  
  /** @brief Scattering cross-section.
   * @param domain Domain over which to compute @f$\sigma_{sc}@f$.
   * @return Scattering cross section and near-field intensity.
   */
  std::vector<double> scatteringCS(Domain* domain);
  
  /** @brief Absorption cross-section.
   * @param domain Domain over which to compute @f$\sigma_{\alpha}@f$.
   * @return Absorption cross section.
   * @note Not always accurate (volumetric methods based on Ohmic losses have better accuracy).
   */
  double absorptionCS(Domain* domain);
  
  /** @brief Extinction cross-section.
   * @param domain Domain over which to compute @f$\sigma_{ext}@f$.
   * @return Extinction cross section.
   * @note Not always accurate (volumetric methods based on Ohmic losses have better accuracy).
   */
  double extinctionCS(Domain* domain);

  // ---------------- Field evaluation over many points ----------------
  /**
   * @brief Evaluate fields on a set of points read from a file.
   *
   * Writes one or more output files (depending on @p fields bitmask):
   *  - .ein/.esc, .hin/.hsc, .din/.dsc, .bin/.bsc, .pin (Poynting), .csc (CS).
   *
   * @param jobFile   Job file path.
   * @param posFile   Position list (CSV of x y z) or empty to skip.
   * @param outPrefix Prefix for output filenames.
   * @param fields    Bitmask: EFIELD/HFIELD/DFIELD/BFIELD/POYVEC/POWER.
   * @param needsCS   If true, writes cross sections to .csc.
   * @return 0 on success; non-zero on failure.
   */
  int FieldEval(std::string jobFile, std::string posFile, std::string outPrefix,
                uint fields, bool needsCS = false);

  /**
   * @brief Thread worker to evaluate fields over @p posvec in parallel.
   * @param posvec   Input positions.
   * @param needsE   Whether to compute @f$\mathbf{E}@f$.
   * @param needsH   Whether to compute @f$\mathbf{H}@f$.
   * @param result   Output array of 
   *                 @f$\{\{\mathbf{E}_\mathrm{inc},  \mathbf{H}_\mathrm{inc}\}, 
   *                      \{\mathbf{E}_\mathrm{scat}, \mathbf{H}_\mathrm{scat}\}\}@f$ 
   *                 per point.
   * @param contDom  Pre-computed containing domains for each point.
   * @param nPoint   Shared index (guarded by @p mut).
   * @param mut      Mutex for index progress + periodic status prints.
   */
  void FieldEvalParallel(
      std::vector<rvec>& posvec, bool needsE, bool needsH,
      std::vector<std::pair<std::pair<cvec, cvec>, std::pair<cvec, cvec>>>&
          result,
      std::vector<Domain*>& contDom, int& nPoint, std::mutex& mut);

  // ---------------- Casimir force computation ----------------
  /** 
   * @brief Casimir workflow entry (uses gradient matrices if enabled). 
   * @param jobFile   Job file path.
   * @param outPrefix Prefix for output filenames.
   * @return 0 on success; non-zero on failure.
   */
  int Casimir(std::string jobFile, std::string outPrefix);

  // ---------------- Utilities ----------------
  /** 
   * @brief Load binary lookup table used by periodic Green's function acceleration. 
   * @param fileName Lookup table file name.
   * @return 0 on success; non-zero on failure.
   */
  int LoadTableFromFile(std::string fileName);

  /** 
   * @brief Set thread count used internally by SimJob. 
   * @param t Thread count.
   */
  static void AssignThreads(int t);

  /** 
   * @brief Dump A into a MATLAB friendly text file.
   * @param lseMatrix Matrix A.
   * @param filename Output filename.
   */
  void dumpLSEMatrix(blitz::Array<dcmplx,2>* lseMatrix, const std::string &filename);

  /** 
   * @brief Dump RHS vector into a MATLAB friendly text file. 
   * @param lseVector RHS vector.
   * @param filename Output filename.
   */
  void dumpLSEVector(blitz::Array<dcmplx,1>* lseVector, const std::string &filename);
};

#endif
