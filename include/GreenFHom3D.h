/**
* @file GreenFHom3D.h
* @brief 3D scalar/dyadic Green's function in a homogeneous background.
*
* Concrete implementation of #GreenF for a single, homogeneous, isotropic
* medium. Provides pointwise scalar, gradient, and dyadic evaluations as well
* as RWG‑weighted integrals used in matrix assembly (D/K) and post‑processing
* (delta/kappa). Includes singularity‑subtraction paths and accuracy control
* via Dunavant quadrature orders.
*/

#ifndef GREENFHOM3D_H
#define GREENFHOM3D_H

#include <blitz/array.h>
#include <cmath>
#include <tuple>
#include "GaussQuad.h"
#include "GreenF.h"
#include "globals.h"

/// \ingroup greenFunction
/**
 * @class GreenFHom3D
 * @brief Free‑space (homogeneous) Green's function provider.
 *
 * Implements:
 * - Scalar Green's function \f$G(r,r') = e^{ik|r-r'|}/(4\pi|r-r'|)\f$.
 * - Gradient \f$\nabla^{\prime} G\f$ and full **dyadic** Green's function.
 * - Smoothed kernels used for singularity subtraction when triangles are
 * close or coincident.
 * - Batched assembly helpers over groups of RWG functions living on the same
 * triangles to reduce repeated kernel evaluations.
 */
class GreenFHom3D : public GreenF {
 private:
  /** @brief Gaussian quadrature helper. */
  GaussQuad<GreenFHom3D> gauss;

 public:
  // ---------------------------------------------------------------------------
  // Proximity tests (heuristics to pick smoothed vs. raw kernels)
  // ---------------------------------------------------------------------------
  /** @brief Heuristic test based on @f$ k |\mathbf{r} − \mathbf{r}^{\prime}| @f$.
   *  @param r Observation point.
   *  @param rp Source point.
   *  @return true if the separation is above the threshold; false otherwise.
   */
  bool AboveThreshold(rvec r, rvec rp);

  /** @brief Heuristic test for two triangles using center distance vs. perimeters.
   *  @param T Observation triangle.
   *  @param Tp Source triangle.
   *  @return true if triangles are sufficiently separated; false if near/coincident. 
   */
  bool AboveThreshold(Triangle* T, Triangle* Tp);

  // ---------------------------------------------------------------------------
  // Construction & lifetime
  // ---------------------------------------------------------------------------
  /** @brief Default constructor (initializes with zeroed wavenumber and no lookup table). */
  GreenFHom3D();
  
  /** @brief Construct for a homogeneous, non-magnetic medium.
   *  @param wavelength Wavelength in the medium.
   *  @param epsilonMedium Relative permittivity \f$\varepsilon_r\f$. 
   */
  GreenFHom3D(dcmplx wavelength, dcmplx epsilonMedium);
  
  /** @brief Construct for a homogeneous, magnetic medium.
   *  @param wavelength Wavelength in the medium.
   *  @param epsilonMedium Relative permittivity \f$\varepsilon_r\f$.
   *  @param muMedium Relative permeability \f$\mu_r\f$. 
   */
  GreenFHom3D(dcmplx wavelength, dcmplx epsilonMedium, dcmplx muMedium);

  /** @brief Virtual destructor. */
  virtual ~GreenFHom3D(){};

  // ---------------------------------------------------------------------------
  // Pointwise kernels
  // ---------------------------------------------------------------------------
  /**
   * @brief Scalar Green's function \f$ G(\mathbf r,\mathbf r^{\prime}) \f$.
   * @param r Observation point \f$\mathbf r\f$.
   * @param rp Source point \f$\mathbf r^{\prime}\f$.
   * @return Complex scalar \f$ G(\mathbf r,\mathbf r^{\prime})\f$.
   */
  virtual dcmplx Evaluate(rvec r, rvec rp);

  /**
   * @brief Dyadic Green's function \f$\overline{\overline G}(\mathbf r,\mathbf r^{\prime})\f$.
   * @param r Observation point \f$\mathbf r\f$.
   * @param rp Source point \f$\mathbf r^{\prime}\f$.
   * @return Complex \f$ 3 \times 3 \f$ dyad.
   */
  virtual cdyad EvaluateDyadic(rvec r, rvec rp);
  
  /**
   * @brief Gradient w.r.t. the source point \f$\nabla^{\prime}G(\mathbf r,\mathbf r^{\prime})\f$.
   * @param r Observation point \f$\mathbf r\f$.
   * @param rp Source point \f$\mathbf r^{\prime}\f$.
   * @return Complex 3-vector \f$\nabla^{\prime} G(\mathbf r,\mathbf r^{\prime})\f$.
   */
  virtual cvec Gradient(rvec r, rvec rp);

  /** @brief Evaluate \f$(G,\ \nabla^{\prime}G,\ \overline{\overline G})\f$ together.
   *  @param r Observation point.
   *  @param rp Source point.
   *  @return Tuple \f$(G,\ \nabla^{\prime}G,\ \overline{\overline G})\f$. 
   */
  virtual std::tuple<dcmplx, cvec, cdyad> EvaluateAll(rvec r, rvec rp);
  
  // ---------------------------------------------------------------------------
  // Smoothed variants used for near‑singular quadrature
  // ---------------------------------------------------------------------------
  /** @brief Smoothed scalar kernel used under singularity subtraction.
   *  @param r Observation point.
   *  @param rp Source point.
   *  @return Complex scalar smoothed kernel. 
   */
  virtual dcmplx Smoothed(rvec r, rvec rp);

  /** @brief Half-smoothed scalar kernel for coincident triangles.
   *  @param r Observation point.
   *  @param rp Source point.
   *  @return Complex scalar half-smoothed kernel. 
   */
  virtual dcmplx halfSmoothed(rvec r, rvec rp);

  /** @brief Smoothed gradient kernel.
   *  @param r Observation point.
   *  @param rp Source point.
   *  @return Complex 3-vector smoothed gradient. 
   */
  virtual cvec GradientSmoothed(rvec r, rvec rp);

  // ---------------------------------------------------------------------------
  // Triangle and post‑processing integrals
  // ---------------------------------------------------------------------------
  /**
   * @brief D-operator matrix element 
   * \f$\langle \mathbf{f}(\mathbf{r}),\,\overline{\overline G},\,\mathbf{f}(\mathbf{r}^{\prime})\rangle\f$.
   * @param f Testing RWG function (on observation triangle).
   * @param fp Basis RWG function (on source triangle).
   * @return Complex scalar entry of \f$\mathbf D\f$.
   * @note Not used here (exits).
   */
  virtual dcmplx IntegrateD(RWGFun* f, RWGFun* fp);

  /**
   * @brief K-operator matrix element 
   * \f$\langle \mathbf{f}(\mathbf{r}),\,\nabla \times \overline{\overline G},\,\mathbf{f}(\mathbf{r}^{\prime})\rangle\f$.
   * @param f Testing RWG function (on observation triangle).
   * @param fp Basis RWG function (on source triangle).
   * @return Complex scalar entry of \f$\mathbf K\f$.
   * @note Not used here (exits).
   */
  virtual dcmplx IntegrateK(RWGFun* f, RWGFun* fp);
  
  /**
   * @brief @f$\boldsymbol \delta@f$ term at an observation point.
   * @param r Observation point \f$\mathbf r\f$.
   * @param fp Source RWG function.
   * @return Complex 3-vector contribution \f$\boldsymbol \delta(\mathbf r)\f$.
   */
  virtual cvec IntegrateDelta(rvec r, RWGFun* fp);

  /**
   * @brief @f$\boldsymbol \kappa@f$ term at an observation point.
   * @param r Observation point \f$\mathbf r\f$.
   * @param fp Source RWG function.
   * @return Complex 3-vector contribution \f$\boldsymbol \kappa(\mathbf r)\f$.
   */
  virtual cvec IntegrateKappa(rvec r, RWGFun* fp);
  
  /**
   * @brief Pair of D elements with swapped triangle roles.
   * @param f Testing RWG function.
   * @param fp Basis RWG function.
   * @return \f$( D(f,fp),\, D(fp,f) )\f$.
   */
  virtual std::pair<dcmplx, dcmplx> PairD(RWGFun* f, RWGFun* fp);

  /**
   * @brief Pair of K elements with swapped triangle roles.
   * @param f Testing RWG function.
   * @param fp Basis RWG function.
   * @return \f$( K(f,fp),\, K(fp,f) )\f$.
   */
  virtual std::pair<dcmplx, dcmplx> PairK(RWGFun* f, RWGFun* fp);

  // ---------------------------------------------------------------------------
  // Assembly helpers
  // ---------------------------------------------------------------------------
  /**
   * @brief Compute \f$\mathbf D\f$ and \f$\mathbf K\f$ for all RWG pairs on two triangles.
   * @param fvec RWGs living on the observation triangle.
   * @param fpvec RWGs living on the source triangle.
   * @param DK Output matrix of pairs \f$(D,K)\f$ indexed by \c [i][j].
   * @details Adapts quadrature and analytic corrections depending on triangle
   * adjacency: coincident, edge‑adjacent, vertex‑adjacent, or disjoint.
   */
  virtual void SameTriDK(
      const std::vector<RWGFun*>& fvec, const std::vector<RWGFun*>& fpvec,
      std::vector<std::vector<std::pair<dcmplx, dcmplx>>>& DK);

  /**
   * @brief Same as @ref SameTriDK plus gradients (if provided by the GF).
   * @param fvec RWGs on the observation triangle.
   * @param fpvec RWGs on the source triangle.
   * @param DK Output matrix of \f$(D,K)\f$ pairs.
   * @param DKG Output matrix of gradient pairs per entry.
   * @note Intended for different domains (no singularity subtraction here).
   */
  virtual void SameTriDKG(
      const std::vector<RWGFun*>& fvec, const std::vector<RWGFun*>& fpvec,
      std::vector<std::vector<std::pair<dcmplx, dcmplx>>>& DK,
      std::vector<std::vector<std::pair<cvec, cvec>>>& DKG);

  /**
   * @brief Compute \f$(\boldsymbol \delta,\boldsymbol \kappa)\f$ for all RWGs on a source triangle.
   * @param r Observation point.
   * @param fpvec RWGs living on the source triangle.
   * @param DeltaKappa Output vector of pairs \f$(\boldsymbol \delta,\boldsymbol \kappa)\f$.
   */
  virtual void SameTriDeltaKappa(
      rvec r, const std::vector<RWGFun*>& fpvec,
      std::vector<std::pair<cvec, cvec>>& DeltaKappa);
};

#endif
