/**
 * @file DomainHom3DPer.cpp
 * @brief Implementation of the homogeneous, periodically replicated 3D domain.
 */

#include "DomainHom3DPer.h"
#include <iostream>
#include "GreenFHom3DPer2D.h"

// -----------------------------------------------------------------------------
// Construction & lifetime
// -----------------------------------------------------------------------------


/**
 * @details If exactly two lattice vectors are provided, a 2D periodic homogeneous Green's
 * function (#GreenFHom3DPer2D) is instantiated. Otherwise, a warning is printed (hook
 * for future 1D/3D periodic Green's function implementations).
 */
DomainHom3DPer::DomainHom3DPer(SurfaceMesh *inMesh, int inDIndex,
                               dcmplx inEpsilon, dcmplx inMu,
                               dcmplx inVacuumWavelength, rvec blochVector,
                               std::vector<rvec> latticeVectors)
    : DomainHom3D(inMesh, inDIndex, inEpsilon, inMu, inVacuumWavelength, 1) {
  if (latticeVectors.size() == 2) {
    grnFun = new GreenFHom3DPer2D(vacuumWavelength, epsilon, mu, blochVector,
                                  latticeVectors[0], latticeVectors[1]);
  } else {
    std::cout << "Warning: number of lattice vectors different than 2"
              << std::endl;
    // NOTE : This is where the call to periodic Green's functions in 1 or 3
    // directions can be inserted
    // In these cases, the size of latticeVectors would be either 1 or 3
  }
}

DomainHom3DPer::~DomainHom3DPer() { delete grnFun; }
