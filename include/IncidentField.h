/**
 * @file IncidentField.h
 * @brief Abstract interface for incident electromagnetic excitations.
 *
 * This header declares the @c IncidentField base class, which defines the
 * common interface for all excitations used by the SIE formulation (i.e.
 * plane wave, dipole, Gaussian beam). Implementations provide pointwise field
 * evaluation (E and H), layer-aware evaluation helpers, and RWG-based
 * excitation integrals used to assemble right-hand sides.
 *
 * The class also exposes a static switch, @p EnableAccurate, that enables a
 * more accurate (and potentially more expensive) numerical path for the
 * excitation integrals, intended for validation or difficult near-field cases.
 */

#ifndef INCIDENTFIELD_H
#define INCIDENTFIELD_H

#include "RWGFun.h"
#include "globals.h"
#include "LayeredMediaUtils.h"
#include "GreenF.h"

/// \ingroup incidentField
/**
 * @class IncidentField
 * @brief Abstract base for incident fields (excitations) used in SIE.
 *
 * Derived classes must implement both pointwise evaluations and the RWG-based
 * integral projections that populate the excitation vectors. When working in a
 * layered background, implementations typically rely on a layered Green's
 * function provided via @ref setGrnFunLayered.
 */
class IncidentField {
protected:
  /** @brief Global toggle: request accurate/expensive integration path. */
  static bool NeedsAccurate;
 public:
  /** @brief Virtual destructor to allow proper cleanup through base pointer. */
  virtual ~IncidentField();
  
  // ---------------------------------------------------------------------------
  // Controls
  // ---------------------------------------------------------------------------
  /**
   * @brief Enable high-accuracy integration for excitations.
   */
  static void EnableAccurate();

  // ---------------------------------------------------------------------------
  // Pointwise field evaluation
  // ---------------------------------------------------------------------------
  /**
   * @brief Evaluate the incident electric field @f$\mathbf{E}^{inc}(\mathbf{r})@f$.
   * @param position Observation point in Cartesian coordinates.
   * @return Complex 3-vector (Ex, Ey, Ez).
   */
  virtual cvec EvaluateE(rvec position) = 0;
  
  /**
   * @brief Evaluate the incident magnetic field @f$\mathbf{H}^{inc}(\mathbf{r})@f$.
   * @param position Observation point in Cartesian coordinates.
   * @return Complex 3-vector (Hx, Hy, Hz).
   */
  virtual cvec EvaluateH(rvec position) = 0;
  
  /**
   * @brief Evaluate @f$(\mathbf{E}^{inc}, \mathbf{H}^{inc})@f$ in layered media.
   * @param position Observation point in Cartesian coordinates.
   * @return A 2-element vector: @f$(\mathbf{E}^{inc}, \mathbf{H}^{inc})@f$, 
   * each a complex 3-vector.
   */
  virtual std::vector<cvec> EvaluateLayeredEandH(rvec position) = 0;
  
  // ---------------------------------------------------------------------------
  // RWG excitation integrals
  // ---------------------------------------------------------------------------
  /**
   * @brief RWG projection of the incident electric field. Performs the scalar 
   * test integral on the triangle supporting the RWG basis function @p f to 
   * assemble the electric-field excitation entry.
   * @param f Pointer to the RWG test function.
   * @return Complex scalar contribution for the excitation vector.
   */
  virtual dcmplx IntegrateE(RWGFun* f) = 0;

  /**
   * @brief RWG projection of the incident magnetic field.
   * @param f Pointer to the RWG test function.
   * @return Complex scalar contribution for the excitation vector.
   */
  virtual dcmplx IntegrateH(RWGFun* f) = 0;

  /**
   * @brief RWG projection of @f$(\mathbf{E}^{inc}, \mathbf{H}^{inc})@f$ in layered media.
   * @param f Pointer to the RWG test function.
   * @return Two complex scalars: @f$(\mathbf{E}^{exc}, \mathbf{H}^{exc})@f$.
   */
  virtual std::vector<dcmplx> IntegrateLayeredEandH(RWGFun* f) = 0;
  
  // ---------------------------------------------------------------------------
  // Type identification
  // ---------------------------------------------------------------------------
  /** @brief True if this excitation is a PlaneWave. 
   *  @return True if PlaneWave illumination exists.
   */
  virtual bool IsPlaneWave() = 0;

  /** @brief True if this excitation is a Dipole. 
   *  @return True if Dipole illumination exists.
   */
  virtual bool IsDipole() = 0;

  /** @brief True if this excitation is a Gaussian beam. 
   *  @return True if Gaussian illumination exists.
   */
  virtual bool IsGaussian() = 0;

  // ---------------------------------------------------------------------------
  // Layered-media configuration
  // ---------------------------------------------------------------------------
  /**
   * @brief Provide layered-medium Green's function and its tabulation grids.
   * @param fGrnFields Pointer to a Green's-function backend (homogeneous or
   * layered). The callee does not take ownership.
   * @param tabGrids Interpolation grids used by the layered backend.
   */
  virtual void setGrnFunLayered(GreenF *fGrnFields, std::vector<Grid> tabGrids) = 0;
};

#endif
