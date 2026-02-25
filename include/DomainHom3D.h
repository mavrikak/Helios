/**
 * @file DomainHom3D.h
 * @brief Homogeneous 3D simulation domain (derived from Domain).
 *
 * This class models a single region filled with homogeneous, isotropic
 * material characterized by complex relative permittivity @f$\varepsilon_r@f$
 * and permeability @f$\mu_r@f$. It provides a homogeneous 3D Green's function
 * and sets up the PMCHWT SIE formulation.
 */

#ifndef DOMAINHOM3D_H
#define DOMAINHOM3D_H

#include "Domain.h"
#include "GreenFHom3D.h"
#include "SurfaceMesh.h"
#include <algorithm>

/// \ingroup domain
/**
 * @class DomainHom3D
 * @brief Concrete domain for homogeneous 3D media.
 *
 * Responsibilities beyond those of #Domain:
 * - Hold constant material parameters (epsilon, mu).
 * - Build a homogeneous Green's function (#GreenFHom3D).
 * - Provide homogeneous Epsilon/Mu queries.
 * - Instantiate a PMCHWT formulation for assembly.
 * - Implement cross‑section post‑processing via triangle quadrature.
 */
class DomainHom3D : public Domain {
 protected:
  dcmplx epsilon; //!< Complex relative permittivity of the domain.
  dcmplx mu;      //!< Complex relative permeability of the domain.
  
  /**
   * @brief Instantiate the homogeneous Green's function for this domain.
   * @return 0 on success; non-zero otherwise.
   */
  virtual int InitGrnFun();

 public:
  // ---------------------------------------------------------------------------
  // Construction & lifetime
  // ---------------------------------------------------------------------------
  /**
   * @brief Construct a homogeneous domain and immediately create its Green's function.
   * @param inMesh Surface mesh containing region @p inDIndex.
   * @param inDIndex Domain region index in @p inMesh.
   * @param inEpsilon Relative permittivity @f$\varepsilon_r@f$.
   * @param inMu Relative permeability @f$\mu_r@f$.
   * @param inVacuumWavelength Free‑space wavelength @f$\lambda_0@f$.
   */
  DomainHom3D(SurfaceMesh *inMesh, int inDIndex, dcmplx inEpsilon, dcmplx inMu,
              dcmplx inVacuumWavelength);
  
  /**
   * @brief Construct a homogeneous domain, optionally skipping Green's function creation.
   * @param inMesh Surface mesh containing region @p inDIndex.
   * @param inDIndex Domain region index in @p inMesh.
   * @param inEpsilon Relative permittivity @f$\varepsilon_r@f$.
   * @param inMu Relative permeability @f$\mu_r@f$.
   * @param inVacuumWavelength Free‑space wavelength @f$\lambda_0@f$.
   * @param parent If @c true, do not call InitGrnFun() here (used by owners/aggregates).
   */
  DomainHom3D(SurfaceMesh *inMesh, int inDIndex, dcmplx inEpsilon, dcmplx inMu,
              dcmplx inVacuumWavelength, bool parent);
  
  /** @brief Destructor (deletes #grnFun). */
  ~DomainHom3D();
  
  // ---------------------------------------------------------------------------
  // Formulation & materials
  // ---------------------------------------------------------------------------
  /**
   * @brief Initialize the PMCHWT formulation object at the given wavelength.
   * @param inVacWavelength Vacuum wavelength \f$\lambda_0\f$ used to configure the formulation.
   * @return 0 on success; non-zero otherwise.
   */
  int InitFormulation(dcmplx inVacWavelength);

  /**
   * @brief Relative permittivity at a spatial location.
   * @param pos Query location.
   * @return Complex relative permittivity \f$\varepsilon_r(\mathbf r)\f$.
   */
  virtual dcmplx Epsilon(rvec pos);
  
  /**
   * @brief Relative permeability at a spatial location.
   * @param pos Query location.
   * @return Complex relative permeability \f$\mu_r(\mathbf r)\f$.
   */
  virtual dcmplx Mu(rvec pos);

  /** @brief Is this a layered domain? Must be defined by subclasses.
   *  @return true for layered domains; false otherwise. 
   */
  virtual bool IsLayered() { return false; }

  // ---------------------------------------------------------------------------
  // Layered‑media helpers
  // ---------------------------------------------------------------------------
  /**
   * @brief Accessor to the Green's function used for layered-media post-processing.
   * @return Non-owning pointer to the layered Green's function instance.
   */
  virtual GreenF* GetGrnFunLayeredFields() const { return nullptr; }
  
  /**
   * @brief Build or refresh layered Green's function tables for positions of interest.
   * @param posvec Vector of observation positions used to seed/refresh tables.
   */
  virtual void newGrnFunLayered(std::vector<rvec> &posvec);

  // ---------------------------------------------------------------------------
  // Cross‑sections (post‑processing)
  // ---------------------------------------------------------------------------
  /**
   * @brief Scattering cross section via surface quadrature.
   * @param solVector Solution coefficients \f$\mathbf{x}\f$.
   * @return Scattering cross section and near-field intesity.
   */
  virtual std::vector<double> scatteringCS(blitz::Array<dcmplx, 1> *solVector);
  
  /**
   * @brief Absorption cross section via surface quadrature.
   * @param solVector Solution coefficients \f$\mathbf{x}\f$.
   * @return Absorption cross section.
   * @note Not always accurate (volumetric methods based on Ohmic losses have better accuracy).
   */
  virtual double absorptionCS(blitz::Array<dcmplx, 1> *solVector);
  
  /**
   * @brief Extinction cross section via surface quadrature.
   * @param solVector Solution coefficients \f$\mathbf{x}\f$.
   * @return Extinction cross section.
   * @note Not always accurate (volumetric methods based on Ohmic losses have better accuracy).
   */
  virtual double extinctionCS(blitz::Array<dcmplx, 1> *solVector);

  // ---------------------------------------------------------------------------
  // Incident field factories
  // ---------------------------------------------------------------------------
  /**
   * @brief Add a PlaneWave excitation (ownership transferred to Domain).
   * @param wavelength Vacuum wavelength \f$\lambda_0\f$ (complex allowed).
   * @param propagationDirection Unit vector \f$\hat{\mathbf k}\f$ (direction of travel).
   * @param polarization Complex polarization vector (transverse to \f$\hat{\mathbf k}\f$).
   * @return 0 on success; non-zero otherwise.
   */
  virtual int AddPlaneWave(dcmplx wavelength, rvec propagationDirection,
  cvec polarization);
  
  /**
   * @brief Add a Dipole excitation (ownership transferred to Domain).
   * @param position Dipole position.
   * @param polarization Complex dipole moment vector.
   * @return 0 on success; non-zero otherwise.
   */
  virtual int AddDipole(rvec position, cvec polarization);
  
  /**
   * @brief Add a Gaussian beam excitation.
   * @param wavelength Vacuum wavelength.
   * @param waist Beam waist radius \f$ w_0 \f$.
   * @param focal_point Beam focus position.
   * @param propagationdirection Unit vector of beam propagation.
   * @param polarization Complex polarization vector.
   * @return 0 on success; non-zero otherwise.
   */
  virtual int AddGaussian(double wavelength, double waist, rvec focal_point,
  rvec propagationdirection, cvec polarization);
};

#endif
