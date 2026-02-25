/**
 * @file SimJob.cpp
 * @brief Implementation for the SimJob class.
 */

#include "SimJob.h"
#include <iostream>
#include <map>
#include <mutex>
#include <thread>
#include "iofunctions.h"

// Be sure to include all possible domain types here:
#include "DomainHom3D.h"
#include "DomainLayered3D.h"
#include "DomainHom3DPer.h"

// Incident types go here:
#include "IncidentField.h"
#include "PlaneWave.h"

// Construct an empty job (no mesh or domains yet).
SimJob::SimJob() : isCasimirJob(false), thereIsALookupTable(false) {}

// Destructor: frees owned domains and lookup table if present.
SimJob::~SimJob() {
  std::vector<Domain *>::iterator domIter;
  for (domIter = domains.begin(); domIter != domains.end(); domIter++)
    delete *domIter;
  delete t;
}

/**
 * @details Parse a full @f$<job>...</job>@f$ block (mesh, domains, wavelength, etc.).
 *          Handles both periodic and non-periodic jobs and ensures consistency
 *          (e.g., number of domains <= mesh->MaxDomainIndex()).
 */
int SimJob::LoadFromFile() {
  bool hasMesh(false);
  bool hasLabel(false);
  jobParser.readTag();
  int domainCounter(0);
  if (jobParser.getTag() == "job") {
    std::string symStatus;
    while (jobParser.lastRead() != "/job") {
      jobParser.readTag();
      if (jobParser.getTag() == "periodic") {
        if (isCasimirJob) {
          std::cout << "Error: Casimir force computation not implemented for "
                       "periodic structures."
                    << std::endl;
          return 1;
        }
        symStatus = jobParser.getTag();
        LoadPeriodicity();
      }

      if (jobParser.getTag() == "label") {
        hasLabel = !LoadLabel();
      }
      if (jobParser.getTag() == "mesh") {
        if (symStatus == "periodic") {
          hasMesh = !LoadMeshPer();
        } else {
          hasMesh = !LoadMesh();
        }
        if (!hasMesh) {
          std::cout << "Error: job missing mesh." << std::endl;
          return 1;
        }
      }
      if (jobParser.getTag() == "wavelength") {
        LoadWavelength();
      }
      if (jobParser.getTag() == "domain") {
        ++domainCounter;
        if (domainCounter > mesh->MaxDomainIndex()) {
          std::cout << "Error: too many domains in job file." << std::endl;
          return 1;
        }
        if (symStatus == "periodic") {
          LoadDomainPer();
        } else {
          LoadDomain();
        }
      }
      if (jobParser.getTag() == "domainLayers") {
        ++domainCounter;
        if (domainCounter > mesh->MaxDomainIndex()) {
          std::cout << "Error: too many domains in job file." << std::endl;
          return 1;
        }
        LoadDomainLayered();
      }
      if (jobParser.getTag() == "periodicdomain") {
        ++domainCounter;
        if (domainCounter > mesh->MaxDomainIndex()) {
          std::cout << "Error: too many domains in job file." << std::endl;
          return 1;
        }
        if (symStatus == "periodic") {
          LoadDomainPer();
        } else {
          std::cout << "Error: periodic domain in non-periodic job."
                    << std::endl;
          return 1;
        }
      }
      if (jobParser.getTag() == "nonperiodicdomain") {
        ++domainCounter;
        if (domainCounter > mesh->MaxDomainIndex()) {
          std::cout << "Error: too many domains in job file." << std::endl;
          return 1;
        }
        LoadDomain();
      }
      if (jobParser.getTag() == "planewave") {
        if (isCasimirJob) {
          std::cout << "Error: Casimir force computation does not support "
                       "illumination."
                    << std::endl;
          return 1;
        }
        LoadPlanewave();
      }
      if (jobParser.getTag() == "dipole") {
        if (isCasimirJob) {
          std::cout << "Error: Casimir force computation does not support "
                       "illumination."
                    << std::endl;
          return 1;
        }
        LoadDipole();
      }
      if (jobParser.getTag() == "gaussian") {
        if (isCasimirJob) {
          std::cout << "Error: Casimir force computation does not support "
                       "illumination."
                    << std::endl;
          return 1;
        }
        LoadGaussian();
      }
    }
  }

  if (!hasLabel) {
    std::cout << "Error: job missing label." << std::endl;
    return 1;
  }
  if (domainCounter < mesh->MaxDomainIndex()) {
    std::cout << "Error: too many domains in mesh file." << std::endl;
    return 1;
  }
  std::cout << "Finished reading job" << std::endl;
  return 0;
}

// Parse one @f$<subjob>...</subjob>@f$ section (label + incident field).
int SimJob::LoadSubjobFromFile() {
  bool hasLabel(false);
  jobParser.readTag();
  if (jobParser.getTag() == "subjob") {
    std::cout << "Reading subjob ..." << std::endl;
    while (jobParser.lastRead() != "/subjob") {
      jobParser.readTag();
      if (jobParser.getTag() == "label") {
        hasLabel = !LoadLabel();
      }
      if (jobParser.getTag() == "planewave") {
        LoadPlanewave();
      }
      if (jobParser.getTag() == "dipole") {
        LoadDipole();
      }
      if (jobParser.getTag() == "gaussian") {
        LoadGaussian();
      }
    }
    if (!hasLabel) {
      std::cout << "Subjob missing label!" << std::endl;
      return 1;
    }
    std::cout << "Finished reading subjob." << std::endl;
  } else
    return 1;
  return 0;
}

// Read @f$<label>@f$ and push to labels vector.
int SimJob::LoadLabel() {
  labels.push_back(jobParser.read<std::string>());
  std::cout << "Job name: " << labels.back() << std::endl;
  jobParser.readTag();
  return 0;
}

// Load a surface mesh (non-periodic).
int SimJob::LoadMesh() {
  mesh = new SurfaceMesh;
  int err(mesh->LoadFromFile(jobParser.read<std::string>()));
  if (err == 1) return err;
  std::cout << "Expecting " << mesh->MaxDomainIndex() << " domains."
            << std::endl;
  jobParser.readTag();
  return 0;
}

// Load a surface mesh and register periodic translations.
int SimJob::LoadMeshPer() {
  mesh = new SurfaceMesh;
  for (unsigned int i(0); i < latticeVectors.size(); ++i) {
    mesh->NewTranslation(latticeVectors[i]);
    mesh->NewTranslation(-latticeVectors[i]);
  }
  int err(mesh->LoadFromFile(jobParser.read<std::string>()));
  if (err == 1) return err;
  std::cout << "Expecting " << mesh->MaxDomainIndex() << " domains."
            << std::endl;
  jobParser.readTag();
  return 0;
}

// Read vacuum wavelength λ0.
int SimJob::LoadWavelength() {
  vacuumWavelength = jobParser.read<dcmplx>();
  std::cout << "Vacuum wavelength: " << vacuumWavelength << "nm." << std::endl;
  jobParser.readTag();
  return 0;
}

// Add a homogeneous domain (ε, μ).
int SimJob::LoadDomain() {
  dcmplx epsilon(jobParser.read<dcmplx>());
  dcmplx mu(jobParser.read<dcmplx>());
  Domain *ptr;
  std::cout << "Domain " << domains.size() << ": dielectric permittivity "
            << epsilon << " and magnetic permeability " << mu << std::endl;
  ptr = new DomainHom3D(mesh, domains.size(), epsilon, mu, vacuumWavelength);
  domains.push_back(ptr);
  jobParser.readTag();
  return 0;
}

// Add a layered domain (@f$<domainLayers>...<layer>...@f$ entries).
int SimJob::LoadDomainLayered() {
  jobParser.readTag();
  std::vector<dcmplx> epsilonIn;
  std::vector<dcmplx> muIn;
  std::vector<double> zValsIn;
  while (jobParser.lastRead() != "/domainLayers") {
    if (jobParser.getTag() == "layer") {
      dcmplx epsilon(jobParser.read<dcmplx>());
      dcmplx mu(jobParser.read<dcmplx>());
      epsilonIn.push_back(epsilon);
      muIn.push_back(mu);
      jobParser.readTag();
    }
    if (jobParser.getTag() == "interface") {
      double zValIn(jobParser.read<double>());
      zValsIn.push_back(zValIn);
      jobParser.readTag();
    }
    jobParser.readTag();
  }
  Domain *ptr;
  std::cout << "Layered domain (" << domains.size() << ") bottom to top layers :" << std::endl;
  for (int i = 0; i < epsilonIn.size(); i++)
  {
    std::cout << "->Layer " << i + 1 << ":"
              << " dielectric permittivity = " << epsilonIn.at(i) 
              << " and magnetic permeability = " << muIn.at(i) << std::endl;
    if (i < epsilonIn.size() - 1) {
      std::cout << "--Interface: z = " << zValsIn.at(i) << std::endl;
    }
  }
  ptr = new DomainLayered3D(mesh, domains.size(), epsilonIn, muIn, zValsIn, vacuumWavelength);
  domains.push_back(ptr);
  jobParser.readTag();
  return 0;
}

// Add a periodic domain block (consistency-checked).
int SimJob::LoadDomainPer() {
  dcmplx epsilon(jobParser.read<dcmplx>());
  dcmplx mu(jobParser.read<dcmplx>());
  Domain *ptr;
  std::cout << "Domain " << domains.size() << ": dielectric permittivity "
            << epsilon << " and magnetic permeability " << mu << std::endl;
  if (latticeVectors.size() != 0) {
    ptr = new DomainHom3DPer(mesh, domains.size(), epsilon, mu,
                             vacuumWavelength, blochVector, latticeVectors);
  } else {
    std::cout << "Error in lattice definition" << std::endl;
    return 1;
  }

  if (thereIsALookupTable) {
    ptr->InitLookupTable(t);
  }
  domains.push_back(ptr);
  jobParser.readTag();
  return 0;
}

// Read @f$<planewave>@f$ block and attach to a domain.
int SimJob::LoadPlanewave() {
  jobParser.readTag();
  rvec propagation(0, 0, 0);
  cvec polarization(0, 0, 0);
  int domainIndex(0);
  while (jobParser.lastRead() != "/planewave") {
    if (jobParser.getTag() == "propagation") {
      double f1(jobParser.read<double>());
      double f2(jobParser.read<double>());
      double f3(jobParser.read<double>());
      propagation = f1, f2, f3;
      propagation /= sqrt(dot(propagation, propagation));
      jobParser.readTag();
    }
    if (jobParser.getTag() == "polarization") {
      dcmplx c1(jobParser.read<dcmplx>());
      dcmplx c2(jobParser.read<dcmplx>());
      dcmplx c3(jobParser.read<dcmplx>());
      polarization = c1, c2, c3;
      jobParser.readTag();
    }
    if (jobParser.getTag() == "domain") {
      int din(jobParser.read<int>());
      domainIndex = din;
      jobParser.readTag();
    }
    jobParser.readTag();
  }
  domains[domainIndex]->AddPlaneWave(vacuumWavelength, propagation,
                                     polarization);
  std::cout << "Adding plane wave to domain " << domainIndex << std::endl;
  return 0;
}

// Read @f$<dipole>@f$ block and attach to a domain.
int SimJob::LoadDipole() {
  jobParser.readTag();
  rvec location(0, 0, 0);
  cvec polarization(0, 0, 0);
  while (jobParser.lastRead() != "/dipole") {
    if (jobParser.getTag() == "location") {
      double f1(jobParser.read<double>());
      double f2(jobParser.read<double>());
      double f3(jobParser.read<double>());
      location = f1, f2, f3;
      jobParser.readTag();
    }
    if (jobParser.getTag() == "polarization") {
      dcmplx c1(jobParser.read<dcmplx>());
      dcmplx c2(jobParser.read<dcmplx>());
      dcmplx c3(jobParser.read<dcmplx>());
      polarization = c1, c2, c3;
      jobParser.readTag();
    }
    jobParser.readTag();
  }
  Domain *inDom(ContainingDomain(location));
  inDom->AddDipole(location, polarization);
  int domIndex(inDom->Index());
  std::cout << "Adding dipole to domain " << domIndex << std::endl;
  return 0;
}

// Read @f$<gaussian>@f$ block and attach to a domain.
int SimJob::LoadGaussian() {
  jobParser.readTag();
  rvec propagation{0., 0., 1.};
  cvec polarisation{1., 0., 0.};
  rvec focus_position{0., 0., 0.};
  double waist = real(vacuumWavelength);
  int domainIndex(0);
  while (jobParser.lastRead() != "/gaussian") {
    if (jobParser.getTag() == "propagation") {
      double p1(jobParser.read<double>());
      double p2(jobParser.read<double>());
      double p3(jobParser.read<double>());
      propagation = p1, p2, p3;
      jobParser.readTag();
    }

    if (jobParser.getTag() == "polarization") {
      dcmplx p1(jobParser.read<dcmplx>());
      dcmplx p2(jobParser.read<dcmplx>());
      dcmplx p3(jobParser.read<dcmplx>());
      polarisation = p1, p2, p3;
      jobParser.readTag();
    }

    if (jobParser.getTag() == "focus_position") {
      double p1(jobParser.read<double>());
      double p2(jobParser.read<double>());
      double p3(jobParser.read<double>());
      focus_position = p1, p2, p3;
      jobParser.readTag();
    }

    if (jobParser.getTag() == "waist_size") {
      double w(jobParser.read<double>());
      waist = w;
      jobParser.readTag();
    }

    if (jobParser.getTag() == "domain") {
      int din(jobParser.read<int>());
      domainIndex = din;
      jobParser.readTag();
    }
    jobParser.readTag();
  }

  domains[domainIndex]->AddGaussian(real(vacuumWavelength), waist,
                                    focus_position, propagation, polarisation);
  std::cout << "Adding paraxial Gaussian to the domain " << domainIndex
            << std::endl;
  return 0;
}

// Read @f$<periodic>@f$ (Bloch vector + lattice vectors) for periodic jobs.
int SimJob::LoadPeriodicity() {
  jobParser.readTag();
  while (jobParser.lastRead() != "/periodic") {
    if (jobParser.getTag() == "bloch") {
      double f1(jobParser.read<double>());
      double f2(jobParser.read<double>());
      double f3(jobParser.read<double>());
      blochVector = f1, f2, f3;
      jobParser.readTag();
    }
    if (jobParser.getTag() == "lattice") {
      double f1(jobParser.read<double>());
      double f2(jobParser.read<double>());
      double f3(jobParser.read<double>());
      rvec lat(f1, f2, f3);
      latticeVectors.push_back(lat);
      jobParser.readTag();
    }
    jobParser.readTag();
  }
  return 0;
}

// Allocate and zero RHS vector; append to lseVectors.
int SimJob::CreateLSEVec() {
  blitz::Array<dcmplx, 1> *vector =
      new blitz::Array<dcmplx, 1>(2 * mesh->EdgeCount());
  *vector = 0;
  lseVectors.push_back(vector);
  return 0;
}

// Allocate and zero solution vector; append to solVectors.
int SimJob::CreateSolVec() {
  blitz::Array<dcmplx, 1> *vector =
      new blitz::Array<dcmplx, 1>(2 * mesh->EdgeCount());
  *vector = 0;
  solVectors.push_back(vector);
  return 0;
}

// Allocate Z (and ∂Z/∂p_i if Casimir path); set to zero.
int SimJob::CreateLSEMatrix() {
  lseMatrix =
      new blitz::Array<dcmplx, 2>(2 * mesh->EdgeCount(), 2 * mesh->EdgeCount(),
                                  blitz::ColumnMajorArray<2>());
  *lseMatrix = 0;
  if (isCasimirJob) {
    for (int i = 0; i < 3; ++i) {
      lseGradientMatrix[i] = new blitz::Array<dcmplx, 2>(
          2 * mesh->EdgeCount(), 2 * mesh->EdgeCount(),
          blitz::ColumnMajorArray<2>());
      *(lseGradientMatrix[i]) = 0;
    }
  }
  return 0;
}

// External conjugate gradient and LU routines
extern "C" {
  /**
  * @brief External conjugate gradient solver (complex version) from cg.f.
  *
  * Solves a linear system Ax = y iteratively using the conjugate gradient
  * squared (CGS) algorithm for complex-valued matrices.
  *
  * @param ca Pointer to the complex system matrix A (flattened).
  * @param cy Pointer to the complex right-hand side vector y.
  * @param cx Pointer to the complex solution vector x (input guess, then result).
  * @param mxn Total number of elements in A.
  * @param n Dimension of the linear system.
  * @param err Target relative error tolerance.
  * @param mxiter Maximum number of iterations allowed.
  * @param iter On exit, the number of iterations performed.
  * @param time On exit, elapsed time in seconds.
  */
  void cgsqr_(dcmplx *ca, dcmplx *cy, dcmplx *cx, int *mxn, int *n, double *err,
              int *mxiter, int *iter, float *time);

  /**
  * @brief Computes the LU factorization of a general M-by-N matrix.
  *
  * @param dim1 Number of rows of the matrix (for `zgetrf_`).
  * @param dim2 Number of columns of the matrix (for `zgetrf_`).
  * @param a Pointer to the matrix to be factored or solved.
  * @param lda Leading dimension of the matrix.
  * @param ipiv Pointer to the pivot indices.
  * @param info Output status (0 if successful).
  */
  void zgetrf_(int *dim1, int *dim2, dcmplx *a, int *lda, int *ipiv, int *info);

  /**
  * @brief Solves a system of linear equations using the LU factorization.
  *
  * @param TRANS Specifies the form of the system of equations ('N' for no transpose, etc.).
  * @param N Order of the coefficient matrix.
  * @param NRHS Number of right-hand sides.
  * @param A Pointer to the LU-factored coefficient matrix.
  * @param LDA Leading dimension of A.
  * @param IPIV Pointer to the pivot indices from `zgetrf_`.
  * @param B Pointer to the right-hand side vector(s).
  * @param LDB Leading dimension of B.
  * @param INFO Output status (0 if successful).
  */
  void zgetrs_(char *TRANS, int *N, int *NRHS, dcmplx *A, int *LDA, int *IPIV,
              dcmplx *B, int *LDB, int *INFO);
}

/**
 * @brief Apply Jacobi preconditioning to a complex linear system.
 *
 * @param inMatrix Pointer to the complex system matrix (modified in place).
 * @param inVector Pointer to the complex right-hand side vector (modified in place).
 * @return Always returns 0.
 * @details Scales each row of the input matrix and corresponding vector element by the
 * inverse of the diagonal entry, effectively normalizing the system's diagonal
 * to unity. For each row i, this routine multiplies both the i-th vector element and the
 * entire row of @p inMatrix by 1 / A(i,i).
 */
int JacobiPrecondition(blitz::Array<dcmplx, 2> *inMatrix,
                       blitz::Array<dcmplx, 1> *inVector) {
  int size(inVector->shape()(0) / 2);
  for (int i = 0; i < 2 * size; i++) {
    dcmplx multElem(1. / (*inMatrix)(i, i));
    (*inVector)(i) *= multElem;
    for (int j = 0; j < 2 * size; j++) {
      (*inMatrix)(i, j) *= multElem;
    }
  }
  return 0;
}

/**
 * @details Solve linear system with selected solver backend.
 *        (cgsqr / LU via external FORTRAN or BLAS/LAPACK calls.)
 */
int SimJob::SolveLSE(int solverType) {
  int error(0);
  if (solverType == 1) {
    int n(2 * mesh->EdgeCount()), maxIter(10000), nIter;
    double err(1e-6);
    float time(0);

    std::cout << "Preconditioning...";
    std::cout.flush();
    JacobiPrecondition(lseMatrix, lseVectors.front());

    std::cout << "Solving system with iterative solver...";
    std::cout.flush();
    cgsqr_(&((*lseMatrix)(0, 0)), &((*(lseVectors.front()))(0)),
           &((*(solVectors.front()))(0)), &n, &n, &err, &maxIter, &nIter,
           &time);
    std::cout << "finished in " << nIter << " steps." << std::endl;
    if (nIter > maxIter) error = 1;
  }
  else if (solverType == 0 || solverType == 2) {
    std::cout << "Solving system with LU decomposition..." << std::endl;
    int n(2 * mesh->EdgeCount());
    blitz::Array<int, 1> iPiv(n);
    int info(0);
    int nRSides(lseVectors.size());
    char trans = 'N';
    blitz::Array<dcmplx, 2> rightSides(n, nRSides,
                                       blitz::ColumnMajorArray<2>());
    for (int i(0); i < nRSides; i++) {
      rightSides(blitz::Range(0, n - 1), i) = *(lseVectors.at(i));
    }
    // Use lapack routines
    zgetrf_(&n, &n, &((*lseMatrix)(0, 0)), &n, &(iPiv(0)), &info);
    zgetrs_(&trans, &n, &nRSides, &((*lseMatrix)(0, 0)), &n, &(iPiv(0)),
            &(rightSides(0, 0)), &n, &info);
    for (int i(0); i < nRSides; i++) {
      *(solVectors.at(i)) = rightSides(blitz::Range(0, n - 1), i);
    }
  }
  return error;
};

// Write a Blitz vector to text.
template <class T>
int SimJob::WriteVectorToFile(std::string fileName,
                              blitz::Array<T, 1> *vector) {
  std::ofstream fsol(fileName.c_str());
  if (!fsol.is_open()) {
    std::cout << "Error opening output file: " << fileName << std::endl;
    return 1;
  }
  fsol << std::scientific;
  fsol.precision(20);
  for (int i = 0; i < vector->shape()(0); i++) {
    fsol << (*vector)(i) << std::endl;
  }
  fsol.close();
  return 0;
}

// Load a saved solution (.sol) into solVectors.back().
int SimJob::LoadSolFromFile(std::string solFile) {
  std::ifstream fsol(solFile.c_str());
  if (!fsol.is_open()) {
    std::cout << "Error opening solution file: " << solFile << std::endl;
    return 1;
  }
  for (int i = 0; i < solVectors.back()->shape()(0); i++) {
    dcmplx buf;
    int err = ReadCommented<dcmplx>(&fsol, &buf);
    if (err) {
      std::cout << "Error reading solution file: " << solFile << " at position "
                << i << std::endl;
      fsol.close();
      return 1;
    }
    (*solVectors.front())(i) = buf;
  }
  fsol.close();
  return 0;
}

// Casimir force workflow
int SimJob::Casimir(std::string jobFile, std::string outPrefix) {
  isCasimirJob = true;

  // Load primary job
  jobParser.open(jobFile);
  if (LoadFromFile()) return 1;

  // Initialize formulations and fill LSE matrix
  CreateLSEMatrix();
  std::vector<Domain *>::iterator domIter;
  for (domIter = domains.begin(); domIter != domains.end(); domIter++) {
    (*domIter)->InitFormulation(vacuumWavelength);
    (*domIter)->FillLSEMatrixAndGradient(lseMatrix, lseGradientMatrix);
  }

  // Finished reading jobs
  jobParser.close();

  int n(2 * mesh->EdgeCount());
  blitz::Array<int, 1> iPiv(n);
  int info(0);
  zgetrf_(&n, &n, &((*lseMatrix)(0, 0)), &n, &(iPiv(0)), &info);
  char trans = 'N';
  blitz::Array<dcmplx, 1> force(3);
  blitz::Array<dcmplx, 1> RWGforce(n);
  std::string fileprefix = outPrefix + labels.at(0) + ".f";
  for (int i = 0; i < 3; ++i) {
    zgetrs_(&trans, &n, &n, &((*lseMatrix)(0, 0)), &n, &(iPiv(0)),
            &((*(lseGradientMatrix[i]))(0, 0)), &n, &info);
    for (int j = 0; j < n; ++j) {
      RWGforce(j) = (*(lseGradientMatrix[i]))(j, j) / (2. * PI);
      force(i) += RWGforce(j);
    }
    char c='x'+i;
    WriteVectorToFile(fileprefix + c, &RWGforce);
  }
  WriteVectorToFile(fileprefix + "xyz", &force);

  // Clean up LSE
  delete lseMatrix;
  while (lseVectors.size() > 0) {
    delete lseVectors.back();
    delete solVectors.back();
    lseVectors.pop_back();
    solVectors.pop_back();
  }
  return 0;
}

// High-level run: parse job, build A and b, solve, and write outputs.
int SimJob::Simulate(std::string jobFile, int solverType,
                     std::string outPrefix) {
  // Load primary job
  jobParser.open(jobFile);
  if (LoadFromFile()) return 1;

  // Initialize formulations and fill LSE matrix
  CreateLSEMatrix();
  std::vector<Domain *>::iterator domIter;
  for (domIter = domains.begin(); domIter != domains.end(); domIter++) {
    (*domIter)->InitFormulation(vacuumWavelength);
    (*domIter)->FillLSEMatrix(lseMatrix);
  }
#ifdef DEBUG
  dumpLSEMatrix(lseMatrix, "LSEMatrix.txt");
#endif
  // Fill LSE vectors for main job and subjobs
  bool validJob(true);
  while (validJob) {
    CreateLSEVec();
    CreateSolVec();
    for (domIter = domains.begin(); domIter != domains.end(); domIter++) {
      (*domIter)->FillLSEVector(lseVectors.back());
      (*domIter)->ClearIncidentFields();
    }
#ifdef DEBUG
    dumpLSEVector(lseVectors.back(), "LSEVector.txt");
    for (int i = 0; i < lseVectors.back()->numElements(); i++) {
      std::cout << "LSE vector element " << i << ": "
                << (*lseVectors.back())(i) << std::endl;
    }
#endif
    validJob = !LoadSubjobFromFile();
  }

  // Finished reading jobs
  jobParser.close();

  // Solve
  if (SolveLSE(solverType) != 0) {
    std::cout << "Iterative solver did not converge, using direct solver!"
              << std::endl;
    SolveLSE(0);
  }

  std::cout << "Writing solution file(s)..." << std::endl;

  for (uint i(0); i < solVectors.size(); i++)
    WriteVectorToFile(
        std::string(outPrefix).append(labels.at(i)).append(".sol"),
        solVectors.at(i));

  // Clean up LSE
  delete lseMatrix;
  while (lseVectors.size() > 0) {
    delete lseVectors.back();
    delete solVectors.back();
    lseVectors.pop_back();
    solVectors.pop_back();
  }
  return 0;
}

// Find the containing domain for pos; returns background if none.
Domain *SimJob::ContainingDomain(rvec pos) {
  for (std::vector<Domain *>::iterator domIter = domains.begin();
       domIter != domains.end(); domIter++)
    if ((*domIter)->IsInside(pos)) return *domIter;
  std::cout << "ERROR! Point " << pos
            << " not found in domains, assuming background!" << std::endl;
  return domains.at(0);
}

// {E_inc, E_scat} without pre-known domain.
std::pair<cvec, cvec> SimJob::EField(rvec pos) {
  Domain *inDom(ContainingDomain(pos));
  if (inDom == NULL)
    return std::pair<cvec, cvec>(cvec(0, 0, 0), cvec(0, 0, 0));
  else
    return std::pair<cvec, cvec>(inDom->IncFieldE(pos),
                                 inDom->SecFieldE(solVectors.back(), pos));
}

// {E_inc, E_scat} wit pre-known domain.
std::pair<cvec, cvec> SimJob::EField(rvec pos, Domain *inDom) {
  return std::pair<cvec, cvec>(inDom->IncFieldE(pos),
                               inDom->SecFieldE(solVectors.back(), pos));
}

// {H_inc, H_scat} without pre-known domain.
std::pair<cvec, cvec> SimJob::HField(rvec pos) {
  Domain *inDom(ContainingDomain(pos));
  if (inDom == NULL)
    return std::pair<cvec, cvec>(cvec(0, 0, 0), cvec(0, 0, 0));
  else
    return std::pair<cvec, cvec>(inDom->IncFieldH(pos),
                                 inDom->SecFieldH(solVectors.back(), pos));
}

// {H_inc, H_scat} with pre-known domain.
std::pair<cvec, cvec> SimJob::HField(rvec pos, Domain *inDom) {
  return std::pair<cvec, cvec>(inDom->IncFieldH(pos),
                               inDom->SecFieldH(solVectors.back(), pos));
}

// {{E_inc,H_inc},{E_scat,H_scat}} at @p pos within @p inDom.
std::pair<std::pair<cvec, cvec>, std::pair<cvec, cvec>> SimJob::EandHField(
    rvec pos, Domain *inDom) {
  return make_pair(make_pair(inDom->IncFieldE(pos), inDom->IncFieldH(pos)),
                   inDom->SecFieldEandH(solVectors.back(), pos));
}

// Derived power using latest solution vector.
double SimJob::Power(Domain *domain) {
  return domain->Power(solVectors.back());
}

// Derived scattering cross section using latest solution vector.
std::vector<double> SimJob::scatteringCS(Domain *domain) {
  return domain->scatteringCS(solVectors.back());
}

// Derived absorption cross section using latest solution vector.
double SimJob::absorptionCS(Domain *domain) {
  return domain->absorptionCS(solVectors.back());
}

// Derived extinction cross section using latest solution vector.
double SimJob::extinctionCS(Domain *domain) {
  return domain->extinctionCS(solVectors.back());
}

/**
 * @brief Close and deallocate a set of open output files.
 *
 * @param files Map of filenames (as C-string keys) to open file stream pointers.
 * @return Always returns 0.
 * @details Iterates through a map of file streams, closing and deleting each 
 * associated @c std::ofstream pointer, then clears the container. Each entry's 
 * stream is closed and its pointer freed using @c delete. The input map is then 
 * cleared to remove all entries.
 */
int CleanUpFiles(std::map<const char *, std::ofstream *> files) {
  for (std::map<const char *, std::ofstream *>::iterator itFile = files.begin();
       itFile != files.end(); itFile++) {
    itFile->second->close();
    delete itFile->second;
  }
  files.clear();
  return 0;
}

/**
 * @brief Compute the element-wise complex conjugate of a 3-vector.
 * @param in Input complex vector.
 * @return Complex vector containing conjugated components of @p in.
 */
cvec ConjVec(cvec in) { return cvec(conj(in(0)), conj(in(1)), conj(in(2))); }

/**
 * @brief Extract the real part of each component of a complex 3-vector.
 * @param in Input complex vector.
 * @return Real vector with the real parts of @p in.
 */
rvec RealVec(cvec in) { return rvec(real(in(0)), real(in(1)), real(in(2))); }

// Batch field evaluation over points from file; writes requested outputs.
int SimJob::FieldEval(std::string jobFile, std::string posFile,
                      std::string outPrefix, uint fields, bool needsCS) {
  // Load primary job
  jobParser.open(jobFile);
  if (LoadFromFile()) return 1;

  CreateSolVec();
  bool validJob(true);
  while (validJob) {
    std::string label = labels.back();
    if(LoadSolFromFile(std::string(label).append(".sol"))) {
      return 1;
    }

    // Map of field files
    std::map<const char *, std::ofstream *> files;

    // Add all needed field files
    if (fields & EFIELD) {
      files["fein"] = new std::ofstream(
          std::string(outPrefix).append(label).append(".ein").c_str());
      files["fesc"] = new std::ofstream(
          std::string(outPrefix).append(label).append(".esc").c_str());
    }
    if (fields & HFIELD) {
      files["fhin"] = new std::ofstream(
          std::string(outPrefix).append(label).append(".hin").c_str());
      files["fhsc"] = new std::ofstream(
          std::string(outPrefix).append(label).append(".hsc").c_str());
    }
    if (fields & DFIELD) {
      files["fdin"] = new std::ofstream(
          std::string(outPrefix).append(label).append(".din").c_str());
      files["fdsc"] = new std::ofstream(
          std::string(outPrefix).append(label).append(".dsc").c_str());
    }
    if (fields & BFIELD) {
      files["fbin"] = new std::ofstream(
          std::string(outPrefix).append(label).append(".bin").c_str());
      files["fbsc"] = new std::ofstream(
          std::string(outPrefix).append(label).append(".bsc").c_str());
    }
    if (fields & POYVEC) {
      files["fpin"] = new std::ofstream(
          std::string(outPrefix).append(label).append(".pin").c_str());
      files["fpsc"] = new std::ofstream(
          std::string(outPrefix).append(label).append(".psc").c_str());
      files["fpex"] = new std::ofstream(
          std::string(outPrefix).append(label).append(".pex").c_str());
    }
    if (fields & POWER) {
      files["fpwr"] = new std::ofstream(
          std::string(outPrefix).append(label).append(".pwr").c_str());
    }

    // Check for output file errors
    bool fileError(false);
    for (std::map<const char *, std::ofstream *>::iterator itFile =
             files.begin();
         itFile != files.end(); itFile++)
      fileError |= !(itFile->second->is_open());
    if (fileError) {
      std::cout << "Error opening output file!" << std::endl;
      CleanUpFiles(files);
      return 1;
    }

    // Calculate forces if desired (no need for points file)
    if (fields & POWER && !needsCS) {
      for (std::vector<Domain *>::iterator iDomain = domains.begin();
           iDomain != domains.end(); iDomain++) {
        double power = Power(*iDomain);
        int index = (*iDomain)->Index();
        *files["fpwr"] << index << ": " << power << std::endl;
      }
    }

    // Are "pointwise" fields desired?
    if (fields - POWER != 0 && !needsCS) {
      std::cout << "checking points file" << std::endl;
      // Check for points file
      std::ifstream fpos(posFile.c_str());
      if (!fpos.is_open()) {
        std::cout << "Error opening input file: " << posFile << std::endl;
        CleanUpFiles(files);
        return 1;
      }

      // Determine required fields
      bool needsE((fields & EFIELD) || (fields & DFIELD) || (fields & POYVEC));
      bool needsH((fields & HFIELD) || (fields & BFIELD) || (fields & POYVEC));

      // Compute and save fields
      std::pair<cvec, cvec> efield;
      std::pair<cvec, cvec> hfield;
      std::vector<rvec> posvec;
      int nPoint(0);
      double x(0), y(0), z(0);
      while (!(ReadCommented<double>(&fpos, &x) +
               ReadCommented<double>(&fpos, &y) +
               ReadCommented<double>(&fpos, &z))) {
        posvec.push_back({x, y, z});
      }
      std::vector<std::pair<std::pair<cvec, cvec>, std::pair<cvec, cvec>>>
          result(posvec.size());
      std::vector<Domain *> contDom(posvec.size());
      // If the domain is layered, create a new tabulated Green's function
      // for fields calculation at the given observation points.
      for (int i = 0; i < posvec.size(); ++i) {
        contDom[i] = ContainingDomain(posvec[i]); assert(contDom[i] != nullptr);
        if (contDom[i]->IsLayered() && contDom[i]->GetGrnFunLayeredFields() == nullptr) {
          contDom[i]->newGrnFunLayered(posvec);
        }
      }
      std::thread *th = new std::thread[threads - 1];
      std::mutex mut;
      for (int i = 0; i < threads - 1; ++i) {
        th[i] = std::thread(&SimJob::FieldEvalParallel, this, std::ref(posvec),
                            needsE, needsH, std::ref(result), std::ref(contDom),
                            std::ref(nPoint), std::ref(mut));
      }
      FieldEvalParallel(posvec, needsE, needsH, result, contDom, nPoint, mut);
      for (int i = 0; i < threads - 1; ++i) {
        th[i].join();
      }
      delete[] th;

      std::cout << "Evaluated " << nPoint << " points." << std::endl;
      std::cout << "Writing fields" << std::endl;
      for (int i = 0; i < nPoint; i++) {
        if (needsE) {
          efield = {result[i].first.first, result[i].second.first};
        }
        if (needsH) {
          hfield = {result[i].first.second, result[i].second.second};
        }
        if (fields & EFIELD) {
          *files["fein"] << efield.first(0) << "  " << efield.first(1) << "  "
                         << efield.first(2) << "  " << std::endl;
          *files["fesc"] << efield.second(0) << "  " << efield.second(1) << "  "
                         << efield.second(2) << "  " << std::endl;
        }
        if (fields & HFIELD) {
          *files["fhin"] << hfield.first(0) << "  " << hfield.first(1) << "  "
                         << hfield.first(2) << "  " << std::endl;
          *files["fhsc"] << hfield.second(0) << "  " << hfield.second(1) << "  "
                         << hfield.second(2) << "  " << std::endl;
        }
        if (fields & DFIELD) {
          dcmplx epsilon = contDom[i]->Epsilon(posvec[i]) * EPS0;
          *files["fdin"] << efield.first(0) * epsilon << "  "
                         << efield.first(1) * epsilon << "  "
                         << efield.first(2) * epsilon << "  " << std::endl;
          *files["fdsc"] << efield.second(0) * epsilon << "  "
                         << efield.second(1) * epsilon << "  "
                         << efield.second(2) * epsilon << "  " << std::endl;
        }
        if (fields & BFIELD) {
          dcmplx mu = contDom[i]->Mu(posvec[i]) * MU0;
          *files["fbin"] << hfield.first(0) * mu << "  " << hfield.first(1) * mu
                         << "  " << hfield.first(2) * mu << "  " << std::endl;
          *files["fbsc"] << hfield.second(0) * mu << "  "
                         << hfield.second(1) * mu << "  "
                         << hfield.second(2) * mu << "  " << std::endl;
        }
        if (fields & POYVEC) {
          rvec poyInc(RealVec(cross(efield.first, ConjVec(hfield.first))) / 2.);
          rvec poySca(RealVec(cross(efield.second, ConjVec(hfield.second))) /
                      2.);
          rvec poyExt(RealVec(cross(efield.first, ConjVec(hfield.second)) +
                              cross(efield.second, ConjVec(hfield.first))) /
                      2.);
          *files["fpin"] << poyInc(0) << "  " << poyInc(1) << "  " << poyInc(2)
                         << "  " << std::endl;
          *files["fpsc"] << poySca(0) << "  " << poySca(1) << "  " << poySca(2)
                         << "  " << std::endl;
          *files["fpex"] << poyExt(0) << "  " << poyExt(1) << "  " << poyExt(2)
                         << "  " << std::endl;
        }
      }
      fpos.close();
    }

    // Calculate scattering, absorption and extinction cross sections
    Domain* background = domains.front();
    if ( needsCS && !( background->hasDipoles() ) ){
      files["fcsc"] = new std::ofstream(
          std::string(outPrefix).append(label).append(".csc").c_str());
      std::vector<double> scaCS = scatteringCS(background);
      double absCS = absorptionCS(background);
      double extCS = extinctionCS(background);
      *files["fcsc"] << scaCS[0] << "  " << absCS << "  " << extCS << "  " << scaCS[1] << std::endl;
    }
    CleanUpFiles(files);

    for (std::vector<Domain *>::iterator domIter = domains.begin();
         domIter != domains.end(); domIter++)
      (*domIter)->ClearIncidentFields();
    background->ClearIncidentFields();
    validJob = !LoadSubjobFromFile();
  }

  return 0;
}

// Worker loop: pulls next point index under a mutex and fills results.
void SimJob::FieldEvalParallel(
    std::vector<rvec> &posvec, bool needsE, bool needsH,
    std::vector<std::pair<std::pair<cvec, cvec>, std::pair<cvec, cvec>>>
        &result,
    std::vector<Domain *> &contDom, int &nPoint, std::mutex &mut) {
  std::pair<cvec, cvec> efield;
  std::pair<cvec, cvec> hfield;
  std::pair<std::pair<cvec, cvec>, std::pair<cvec, cvec>> ehfield;
  while (1) {
    int i;
    {
      // Important : do not remove braces, they are there to block-scope mutex
      // lock
      std::lock_guard<std::mutex> lock(mut);
      i = nPoint;
      if (i >= result.size()) break;
      if (++nPoint % 10 == 0) {
        std::cout << ".";
        std::cout.flush();
      }
    }
    if (needsE && needsH) {
      ehfield = EandHField(posvec[i], contDom[i]);
      result[i] = ehfield;
    } else if (needsE) {
      efield = EField(posvec[i], contDom[i]);
      result[i].first.first = efield.first;
      result[i].second.first = efield.second;
    } else if (needsH) {
      hfield = HField(posvec[i], contDom[i]);
      result[i].first.second = hfield.first;
      result[i].second.second = hfield.second;
    }
  }
}

// Load binary lookup table for periodic Green's function helper.
int SimJob::LoadTableFromFile(std::string fileName) {
  std::cout << "Loading lookup table..." << std::endl;
  t = new LookupTableBin(fileName);
  thereIsALookupTable = true;
  return 0;
}

// Global thread setting for SimJob helper routines.
int SimJob::threads = 1;
void SimJob::AssignThreads(int t) { threads = t; }

// Dump A in a MATLAB friendly format.
void SimJob::dumpLSEMatrix(blitz::Array<dcmplx,2>* lseMatrix,
                           const std::string &filename)
{
    if (!lseMatrix) {
        throw std::runtime_error("lseMatrix pointer is null");
    }

    std::ofstream fout(filename);
    if (!fout.is_open()) {
        throw std::runtime_error("Could not open " + filename + " for writing");
    }

    int M = lseMatrix->extent(0);
    int N = lseMatrix->extent(1);

    fout << "A = [\n";
    for (int i = 0; i < M; ++i) {
        fout << "  ";
        for (int j = 0; j < N; ++j) {
            dcmplx z = (*lseMatrix)(i,j);
            // write as real + imag*1i
            fout 
              << std::real(z)
              << (std::imag(z) < 0 ? "" : "+")
              << std::imag(z) << "*1i";
            if (j < N-1) fout << "  ";  // space between entries
        }
        fout << ";\n";  // end of this row
    }
    fout << "];\n";

    fout.close();
}

// Dump RHS in a MATLAB friendly format.
void SimJob::dumpLSEVector(blitz::Array<dcmplx,1>* lseVector,
                           const std::string &filename)
{
    if (!lseVector) {
        throw std::runtime_error("lseVector pointer is null");
    }

    std::ofstream fout(filename);
    if (!fout.is_open()) {
        throw std::runtime_error("Could not open " + filename + " for writing");
    }

    int M = lseVector->extent(0);

    fout << "v = [ ";
    for (int i = 0; i < M; ++i) {
        dcmplx z = (*lseVector)(i);
        fout
            << std::real(z)
            << (std::imag(z) < 0 ? "" : "+")
            << std::imag(z) << "*1i";
        if (i < M-1) fout << "  ";  // space between entries
    }
    fout << " ];\n";

    fout.close();
}
