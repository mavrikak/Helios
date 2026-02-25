/**
 * @file Gaussian.h
 * @brief Incident field: paraxial Gaussian beam in a homogeneous medium.
 *
 * The #Gaussian class implements a Gaussian beam illumination. It evaluates 
 * complex electric and magnetic fields at arbitrary points and provides 
 * RWG‑weighted integrals used in the MoM/SIE RHS.
 *
 * Notes:
 * - The beam is defined by wavelength, waist size @f$ w_0 @f$, a focal point,
 * medium permittivity (via refractive index), propagation direction, and
 * polarization vector.
 * - The implementation assumes a homogeneous background (no layered Green's
 * functions); layered variants throw.
 */

#ifndef GAUSSIAN_H
#define GAUSSIAN_H

#include <iostream>
#include "GaussQuad.h"
#include "IncidentField.h"
#include "globals.h"

/// \ingroup incidentField
/**
* @class Gaussian
* @brief Paraxial Gaussian beam illumination (homogeneous medium).
*
* Field model:
* @f[
* E(r) = \mathbf p\, \frac{w_0}{w(z)} \exp\!
* \left(-\frac{\rho^2}{w^2(z)} + i k [\,Z + \tfrac{1}{2} C(Z)\,\rho^2\,] - i\,\zeta(Z)\right),
* @f]
* with standard definitions for beam radius @f$ w(z) @f$, curvature @f$ C(Z) @f$,
* Gouy phase @f$\zeta(Z)@f$, and cylindrical radius @f$\rho@f$ measured wrt
* the propagation axis. The magnetic field follows from
* @f$\mathbf H = (1/Z_m)\,\hat k \times \mathbf E@f$.
*/
class Gaussian : public IncidentField {
 public:
  // ---------------------------------------------------------------------------
  // Construction
  // ---------------------------------------------------------------------------
  /**
   * @brief Construct a Gaussian beam in a homogeneous medium.
   * @param wavelength Wavelength in vacuum.
   * @param waist Beam waist @f$ w_0 @f$ at the focus.
   * @param focal_point Focus position.
   * @param epsilonMedium Relative permittivity of the medium.
   * @param propagationdirection Unit vector (not necessarily normalized on input).
   * @param polarisation Complex polarization vector (unit‑like magnitude recommended).
   */
  Gaussian(double wavelength, double waist, rvec focal_point,
           double epsilonMedium, rvec propagationdirection, cvec polarisation);

  // ---------------------------------------------------------------------------
  // Pointwise evaluation
  // ---------------------------------------------------------------------------
  /**
   * @brief Evaluate the incident electric field @f$\mathbf{E}^{inc}(\mathbf{r})@f$.
   * @param r Observation point in Cartesian coordinates.
   * @return Complex 3-vector (Ex, Ey, Ez).
   */
  cvec EvaluateE(rvec r);

  /**
   * @brief Evaluate the incident magnetic field @f$\mathbf{H}^{inc}(\mathbf{r})@f$.
   * @param r Observation point in Cartesian coordinates.
   * @return Complex 3-vector (Hx, Hy, Hz).
   */
  cvec EvaluateH(rvec r);
  
  /**
   * @brief Evaluate @f$(\mathbf{E}^{inc}, \mathbf{H}^{inc})@f$ in layered media.
   * @param r Observation point in Cartesian coordinates.
   * @return A 2-element vector: {E_inc, H_inc}, each a complex 3-vector.
   * @note Not implemented yet.
   */
  std::vector<cvec> EvaluateLayeredEandH(rvec r);

  // ---------------------------------------------------------------------------
  // RWG‑weighted integrals (RHS contributions)
  // ---------------------------------------------------------------------------
  /**
   * @brief RWG projection of the incident electric field.
   * @param f Pointer to the RWG basis function.
   * @return Complex scalar contribution for the excitation vector.
   * @details Performs the scalar test integral on the triangle supporting 
   * the RWG basis function @p f to assemble the electric-field excitation entry.
   */
  dcmplx IntegrateE(RWGFun* f);

  /**
   * @brief RWG projection of the incident magnetic field.
   * @param f Pointer to the RWG basis function.
   * @return Complex scalar contribution for the excitation vector.
   * @details Performs the scalar test integral on the triangle supporting 
   * the RWG basis function @p f to assemble the magnetic-field excitation entry.
   */
  dcmplx IntegrateH(RWGFun* f);

  /**
   * @brief RWG projection of @f$(\mathbf{E}^{inc}, \mathbf{H}^{inc})@f$ in layered media.
   * @param f Pointer to the RWG basis function.
   * @return Two complex scalars: {E_exc, H_exc}.
   * @note Not implemented yet (throws for layered media).
   */
  std::vector<dcmplx> IntegrateLayeredEandH(RWGFun* f);

  // ---------------------------------------------------------------------------
  // Accessors
  // ---------------------------------------------------------------------------
  /** @brief Beam waist \f$ w_0 \f$.
   *  @return Waist radius at focus. 
   */
  double Waist() const;

  /** @brief Polarization vector.
   *  @return Complex polarization \f$\mathbf p\f$.
   */
  cvec Polarisation() const;

  /** @brief Propagation vector \f$ k \f$.
   *  @return Real 3-vector in the propagation direction scaled by \f$ k \f$. 
   */
  rvec Propagation() const;

 private:
  // --- Derived/physical parameters ---
  double k;     //!< Wavenumber in medium.
  rvec prop;    //!< Unit propagation direction.
  cvec p;       //!< Complex polarization vector.
  double w0;    //!< Waist size.
  double depth; //!< Focus depth.
  rvec focus;   //!< Focal point position.
  double z;     //!< Medium impedance.

  GaussQuad<Gaussian> gauss;  //!< Triangle quadrature helper used by IntegrateE/H.
  
  // --- IncidentField type tags ---
  /** @brief Type tag: this incident field is not a PlaneWave.
   *  @return false. 
   */
  bool IsPlaneWave() { return false; }

  /** @brief Type tag: this incident field is not a Dipole.
   *  @return false. 
   */
  bool IsDipole()    { return false; }

  /** @brief Type tag: this incident field is a Gaussian beam.
   *  @return true. 
   */
  bool IsGaussian()  { return true;  }
  
  // --- Layered GF wiring (unsupported for this class) ---
  /**
   * @brief Set layered Green's function and tabulation grids.
   * @param fGrnFields Pointer to layered Green's function used for field evaluation.
   * @param tabGrids Precomputed tabulation grids for layered post-processing.
   * @note Not implemented yet (throws for layered media).
   */
  virtual void setGrnFunLayered(GreenF *fGrnFields, std::vector<Grid> tabGrids);
};

#endif
