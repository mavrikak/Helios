/**
 * @file DomainHom3DPer.h
 * @brief Homogeneous 3D domain with in‑plane periodicity (derived from DomainHom3D).
 *
 * This domain models an infinite 2D lattice of identical scatterers (periodic
 * along two lattice vectors lying in a plane). It provides a periodic
 * homogeneous Green's function and reuses the rest of the machinery from
 * #DomainHom3D (RWGs, formulation, etc.). For now, only 2D periodicity is
 * instantiated; 1D/3D hooks can be added similarly.
 */

#ifndef DOMAINHOM3DPER_H
#define DOMAINHOM3DPER_H

#include "DomainHom3D.h"
#include "SurfaceMesh.h"

/// \ingroup domain
/**
 * @class DomainHom3DPer
 * @brief Homogeneous medium with Bloch‑periodic excitation/Green's function.
 *
 * The constructor selects a periodic Green's function according to the number
 * of provided lattice vectors. Currently a 2D periodic Green's function is
 * created when two lattice vectors are given; otherwise a warning is printed
 * (extension point for 1D/3D periodicity).
 */
class DomainHom3DPer : public DomainHom3D {
 public:
  /**
   * @brief Construct a periodic homogeneous domain and select a periodic Green's function.
   *
   * @param inMesh Surface mesh containing the region @p inDIndex.
   * @param inDIndex Region index within @p inMesh.
   * @param inEpsilon Relative permittivity of the homogeneous medium.
   * @param inMu Relative permeability of the homogeneous medium.
   * @param inVacuumWavelength Free‑space wavelength.
   * @param blochVector Bloch vector used for Floquet–Bloch conditions.
   * @param latticeVectors List of lattice vectors defining periodicity. If the list
   * has size 2, a 2D periodic Green's function is created with these vectors.
   */
  DomainHom3DPer(SurfaceMesh *inMesh, int inDIndex, dcmplx inEpsilon,
                 dcmplx inMu, dcmplx inVacuumWavelength,
                 rvec blochVector, std::vector<rvec> latticeVectors);

  /** @brief Destructor: deletes the internally owned Green's function. */
  ~DomainHom3DPer();
};

#endif
