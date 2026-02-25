/**
 * @file Dipole.h
 * @brief Incident field generator for electric dipole illumination in homogeneous or layered media.
 *
 * The **Dipole** class derives from #IncidentField and provides the complex
 * electric and magnetic fields radiated by a point electric dipole. It supports
 * both direct field evaluation at arbitrary observation points and RWG-based
 * integrations needed during MoM/SIE matrix assembly. For layered media, it can
 * use tabulated Sommerfeld integrals via #GreenFLayered3D and helper utilities
 * provided by #LayeredMediaUtils.
 */

#ifndef DIPOLE_H
#define DIPOLE_H

#include "GaussQuad.h"
#include "GreenFLayered3D.h"
#include "IncidentField.h"
#include "globals.h"
#include "LayeredMediaUtils.h"

/// \ingroup incidentField
/**
 * @class Dipole
 * @brief Point electric dipole excitation.
 *
 * The dipole is defined by its complex polarization vector @f$\mathbf p@f$
 * and source position @f$\mathbf r_0@f$.
 *
 * In a homogeneous background, fields are obtained from the dyadic free‑space
 * Green's functions. In a stratified background, reflected/transmitted
 * contributions are added using tabulated layered‑media Green's dyadics.
 */
class Dipole : public IncidentField {
 private:
  // --- Physical dipole parameters ---
  cvec p;            //!< Dipole polarization vector
  rvec r0;           //!< Dipole source location
  dcmplx omega;      //!< Angular frequency
  
  // --- Background medium (homogeneous) ---
  GreenF* grnFun;    //!< Pointer to homogeneous background Green's function provider
  dcmplx muMedium;   //!< Relative magnetic permeability of the background medium

  // --- Tabulation data (layered) ---
  std::vector<Grid> tabulationGrids;       //!< Grids for matrix fill‑in (assembly time)
  std::vector<Grid> tabulationGridsFields; //!< Grids for post‑processing field evaluation
  
  // --- Layered‑medium state ---
  GreenF *grnFunLayered;           //!< Layered Green's function (assembly)
  GreenF *grnFunLayeredFields;     //!< Layered Green's function (field evaluation)
  LayeredMediaUtils* layeredUtils; //!< Utility with layer indices, k_z, Fresnel, etc.

 public:
  // ---------------------------------------------------------------------------
  // Construction
  // ---------------------------------------------------------------------------
  /**
   * @brief Construct a dipole in a homogeneous background.
   * @param location Dipole position @f$\mathbf r_0@f$.
   * @param polarization Complex polarization vector @f$\mathbf p@f$.
   * @param vacuumWavelength Free‑space wavelength @f$\lambda_0@f$.
   * @param muMedium Relative permeability of the background.
   * @param fHom Homogeneous Green function provider used for E/H and integrals.
   *
   * The angular frequency is set to @f$\omega = 2\pi c/\lambda_0@f$.
   */
  Dipole(rvec location, cvec polarization, dcmplx vacuumWavelength,
         dcmplx muMedium, GreenF* fHom);

  /**
   * @brief Construct a dipole prepared for layered‑media calculations.
   * @param location Dipole position @f$\mathbf r_0@f$.
   * @param polarization Complex polarization vector @f$\mathbf p@f$.
   * @param vacuumWavelength Free‑space wavelength @f$\lambda_0@f$.
   * @param fHom Homogeneous Green function.
   * @param muMedium Relative permeability of the local layer.
   * @param tabGrids Precomputed tabulation grids for layered kernels.
   * @param fLayered Layered Green function for assembly.
   * @param utils Layered‑media utilities (layer indices, k_z, Fresnel coeffs, etc.).
   */
  Dipole(rvec location, cvec polarization, dcmplx vacuumWavelength, 
         GreenF* fHom, dcmplx muMedium, std::vector<Grid> tabGrids, 
         GreenF* fLayered, LayeredMediaUtils* utils);

  // ---------------------------------------------------------------------------
  // Pointwise field evaluation
  // ---------------------------------------------------------------------------
  /**
   * @brief Evaluate the incident electric field @f$\mathbf E^{inc}(\mathbf r)@f$ in a homogeneous background.
   * @param position Observation point @f$\mathbf r@f$.
   * @return Complex 3‑vector of the electric field.
   */
  cvec EvaluateE(rvec position);
  
  /**
   * @brief Evaluate the incident magnetic field @f$\mathbf H^{inc}(\mathbf r)@f$ in a homogeneous background.
   * @param position Observation point @f$\mathbf r@f$.
   * @return Complex 3‑vector of the magnetic field.
   */
  cvec EvaluateH(rvec position);

  /**
   * @brief Evaluate @f$(\mathbf E^{inc}, \mathbf H^{inc})@f$ in a stratified background.
   * @param position Observation point @f$\mathbf r@f$.
   * @return Two‑element vector @f$\{ \mathbf{E}, \mathbf{H} \}@f$ with complex 3‑vectors.
   * @throws std::runtime_error if layered Green's function for field evaluation was not set.
   * @note If the observation and source points are in the same layer, the homogeneous
   * contribution is added automatically; layered contributions are added on top.
   */
  std::vector<cvec> EvaluateLayeredEandH(rvec position);
  
  // ---------------------------------------------------------------------------
  // RWG integrations (MoM right‑hand‑sides)
  // ---------------------------------------------------------------------------
  /**
   * @brief Integrate the electric‑field kernel against an RWG basis function (homogeneous).
   * @param f RWG test function pointer (the triangle on which it lives defines the support).
   * @return Scalar complex contribution to the RHS of the final system.
   */
  dcmplx IntegrateE(RWGFun* f);

  /**
   * @brief Integrate the magnetic‑field kernel against an RWG basis function (homogeneous).
   * @param f RWG test function pointer.
   * @return Scalar complex contribution to the RHS of the final system.
   */
  dcmplx IntegrateH(RWGFun* f);

  /**
   * @brief Layered‑media counterparts of #IntegrateE and #IntegrateH.
   * @param f RWG test function pointer.
   * @return Two complex scalars @f$\{E^{inc}, H^{inc}\}@f$ combining direct and layered parts.
   */
  std::vector<dcmplx> IntegrateLayeredEandH(RWGFun* f);
  
  // ---------------------------------------------------------------------------
  // Accessors
  // ---------------------------------------------------------------------------
  /** @brief Get dipole polarization. 
   *  @return Dipole polarization complex vector.
   */
  cvec Polarization();
  
  /** @brief Get dipole position. 
   *  @return Dipole location real vector.
   */
  rvec Location();

  // ---------------------------------------------------------------------------
  // Type traits (for polymorphic #IncidentField handling)
  // ---------------------------------------------------------------------------
  bool IsPlaneWave() { return false; }
  bool IsDipole() { return true; }
  bool IsGaussian() { return false; }

  // ---------------------------------------------------------------------------
  // Layered‑media configuration
  // ---------------------------------------------------------------------------
  /**
   * @brief Provide the layered‑media Green's function for post‑processing fields.
   * @param fGrnFields #GreenF instance used to evaluate layered dyadics at points.
   * @param tabGrids Grids used by @p fGrnFields for interpolation.
   * @note This does not modify the assembly‑time (#grnFunLayered) object.
   */
  virtual void setGrnFunLayered(GreenF *fGrnFields, std::vector<Grid> tabGrids);
};

#endif
