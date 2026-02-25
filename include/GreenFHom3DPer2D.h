/**
 * @file GreenFHom3DPer2D.h
 * @brief Periodic homogeneous 3D Green's function with **2D in‑plane periodicity**.
 *
 * Implements a Bloch‑periodic free‑space Green's function for structures
 * periodic along two lattice vectors (a1, a2) in an otherwise homogeneous
 * medium. Convergence is accelerated via **Ewald splitting** into a spatial
 * (real‑space images) sum and a spectral (reciprocal‑lattice) sum. Both the
 * scalar kernel and its source‑gradient are provided, together with RWG‑weighted
 * triangle integrals used in the SIE formulation and post‑processing.
 */

#ifndef GREENFHOM3DPER2D_H
#define GREENFHOM3DPER2D_H

#include <cmath>
#include "GreenFHom3DPer.h"
#include "globals.h"
#include "mathFunctions.h"

/// \ingroup greenFunction
/**
 * @class GreenFHom3DPer2D
 * @brief 2D‑periodic homogeneous 3D Green's function with Ewald acceleration.
 *
 * Key responsibilities:
 * - Manage lattice/reciprocal vectors and Bloch vector projection into the
 * periodic plane (`ToPlane`).
 * - Perform Ewald‑split evaluations of the scalar kernel `G` and its
 * source‑gradient `grad'G` (both smoothed and unsmoothed variants).
 * - Provide RWG‑weighted integrals for RHS/post processing (delta, kappa) and
 * same‑triangle batched assembly helpers.
 */
class GreenFHom3DPer2D : public GreenFHom3DPer {
 private:
  // ---------------------------------------------------------------------------
  // Quadrature helper
  // ---------------------------------------------------------------------------
  GaussQuad<GreenFHom3DPer2D> gauss; //!< Gaussian quadrature for numerical integration.

  // ---------------------------------------------------------------------------
  // Initialization
  // ---------------------------------------------------------------------------
  /** @brief Common initialization: set lookup flags, build reciprocal lattice, 
   * project k, choose Ewald parameter.
   */
  void Init();

 protected:
  // ---------------------------------------------------------------------------
  // Lattice/reciprocal geometry
  // ---------------------------------------------------------------------------
  rvec a1;              //!< First lattice vector (periodic direction).
  double a1a1;          //!< a1·a1.
  double norma1;        //!< |a1|.
  rvec b1;              //!< First reciprocal lattice vector.

  rvec a2;              //!< Second lattice vector (periodic direction).
  double a2a2;          //!< a2·a2.
  double norma2;        //!< |a2|.
  double a1a2;          //!< a1·a2.
  double crossFactor;   //!< |a1xa2|^2 = a1a1*a2a2 − (a1·a2)^2 (used repeatedly).
  rvec b2;              //!< Second reciprocal lattice vector.
  rvec perp;            //!< Unit vector normal to the lattice plane.
  
  // ---------------------------------------------------------------------------
  // Lattice helpers
  // ---------------------------------------------------------------------------
  /**
   * @brief Map a displacement to the primitive cell.
   * @param target Displacement vector.
   * @return The lattice translation \f$\mathbf t\f$ that was applied to @p target.
   */
  rvec PrimitiveCell(rvec& target);
  
  /**
   * @brief Build reciprocal lattice vectors.
   * @details Computes @c b1 and @c b2 and updates cached geometric factors.
   */
  void Reciprocal();
  
  /**
   * @brief Project the Bloch vector onto the lattice plane.
   * @details Removes the component normal to the plane spanned by @c a1 and @c a2.
   */
  void ToPlane();

  // ---------------------------------------------------------------------------
  // Ewald's acceleration technique
  // ---------------------------------------------------------------------------
  /**
   * @brief Ewald-summed scalar kernel (2D periodic).
   * @param r Displacement vector.
   * @return Complex scalar.
   */
  dcmplx Ewald(rvec r);

  /**
   * @brief Ewald-summed smoothed scalar kernel.
   * @param r Displacement vector.
   * @return Complex scalar (used near singularity).
   */
  dcmplx EwaldSmoothed(rvec r);

  /**
   * @brief Ewald-summed source-gradient.
   * @param r Displacement vector.
   * @return Complex 3-vector gradient.
   */
  cvec GradientEwald(rvec r);

  /**
   * @brief Ewald-summed smoothed source-gradient.
   * @param r Displacement vector.
   * @return Complex 3-vector gradient of the smoothed kernel.
   */
  cvec GradientEwaldSmoothed(rvec r);

  // ---------------------------------------------------------------------------
  // Ewald split — spatial (real‑space images)
  // ---------------------------------------------------------------------------
  /**
   * @brief Spatial (real-space images) part of Ewald sum.
   * @param r Displacement vector.
   * @return Complex scalar spatial contribution.
   */
  dcmplx EwaldSpatial(rvec r);

  /**
   * @brief Spatial part of the source-gradient.
   * @param r Displacement vector.
   * @return Complex 3-vector spatial contribution.
   */
  cvec GradientEwaldSpatial(rvec r);

  /**
   * @brief Spatial part returning both normal and gradient.
   * @param r Displacement vector.
   * @return Pair \f$(G,\ \nabla^{\prime}G)\f$ (spatial contribution only).
   */
  std::pair<dcmplx, cvec> EwaldAndGradientSpatial(rvec r);

  /**
   * @brief Spatial part for the smoothed scalar kernel \f$ G_s \f$.
   * @param r Displacement vector.
   * @return Complex scalar spatial contribution (smoothed).
   */
  dcmplx EwaldSpatialSmoothed(rvec r);

  /**
   * @brief Spatial part of the smoothed source-gradient \f$\nabla^{\prime}G_s\f$.
   * @param r Displacement vector.
   * @return Complex 3-vector spatial contribution (smoothed).
   */
  cvec GradientEwaldSpatialSmoothed(rvec r);

  /**
   * @brief Spatial part returning \f$\{G_s,\ \nabla^{\prime} G_s\}\f$.
   * @param r Displacement vector.
   * @return Pair \f$(G_s,\ \nabla^{\prime} G_s)\f$ (spatial contribution only).
   */
  std::pair<dcmplx, cvec> EwaldAndGradientSpatialSmoothed(rvec r);

  // ---------------------------------------------------------------------------
  // Ewald split — spectral (reciprocal‑lattice sum)
  // ---------------------------------------------------------------------------
  /**
   * @brief Spectral (reciprocal-lattice) part of Ewald sum for \f$ G \f$.
   * @param r Displacement vector.
   * @return Complex scalar spectral contribution.
   */
  dcmplx EwaldSpectral(rvec r);

  /**
   * @brief Spectral part of the source-gradient \f$\nabla^{\prime} G\f$.
   * @param r Displacement vector.
   * @return Complex 3-vector spectral contribution.
   */
  cvec GradientEwaldSpectral(rvec r);

  /**
   * @brief Spectral part returning \f$\{G,\ \nabla^{\prime}G\}\f$.
   * @param r Displacement vector.
   * @return Pair \f$(G,\ \nabla^{\prime}G)\f$ (spectral contribution only).
   */
  std::pair<dcmplx, cvec> EwaldAndGradientSpectral(rvec r);

 public:
  // ---------------------------------------------------------------------------
  // Construction
  // ---------------------------------------------------------------------------
  /**
   * @brief Constructor for a non-magnetic medium (2D periodic).
   * @param wavelength Vacuum wavelength \f$\lambda_0\f$.
   * @param epsilonMedium Relative permittivity \f$\varepsilon_r\f$.
   * @param wavevector Bloch vector \f$\mathbf k_B\f$.
   * @param latticeVector First lattice vector \f$\mathbf a_1\f$.
   * @param latticeVector2 Second lattice vector \f$\mathbf a_2\f$.
   */
  GreenFHom3DPer2D(dcmplx wavelength, dcmplx epsilonMedium, rvec wavevector,
                   rvec latticeVector, rvec latticeVector2);

  /**
   * @brief Constructor for a magnetic medium (2D periodic).
   * @param wavelength Vacuum wavelength \f$\lambda_0\f$.
   * @param epsilonMedium Relative permittivity \f$\varepsilon_r\f$.
   * @param muMedium Relative permeability \f$\mu_r\f$.
   * @param wavevector Bloch vector \f$\mathbf k_B\f$.
   * @param latticeVector First lattice vector \f$\mathbf a_1\f$.
   * @param latticeVector2 Second lattice vector \f$\mathbf a_2\f$.
   */
  GreenFHom3DPer2D(dcmplx wavelength, dcmplx epsilonMedium, dcmplx muMedium,
                   rvec wavevector, rvec latticeVector, rvec latticeVector2);
  
  // ---------------------------------------------------------------------------
  // Thresholds & translations
  // ---------------------------------------------------------------------------
  /**
   * @brief Far/near heuristic for translated triangle pairs.
   * @param T Observation triangle.
   * @param Tp Source triangle.
   * @param t Lattice translation applied to the source.
   * @return true if sufficiently separated; false if near.
   */
  bool AboveThresholdTranslate(Triangle* T, Triangle* Tp, rvec t);
  
  // ---------------------------------------------------------------------------
  // Pointwise evaluations (Bloch‑periodic)
  // ---------------------------------------------------------------------------
  /**
   * @brief Scalar Bloch-periodic Green's function \f$ G(\mathbf r,\mathbf r')\f$.
   * @param r Observation point.
   * @param rp Source point.
   * @return Complex scalar \f$ G \f$.
   */
  dcmplx Evaluate(rvec r, rvec rp);

  /**
   * @brief Scalar \f$ G \f$ for a specific lattice translation.
   * @param r Observation point.
   * @param rp Source point.
   * @param t Lattice translation applied to the source (no primitive-cell mapping).
   * @return Complex scalar \f$ G \f$.
   */
  dcmplx EvaluateTranslate(rvec r, rvec rp, rvec t);

  /**
   * @brief Smoothed scalar kernel \f$ G_s(\mathbf r,\mathbf r^{\prime})\f$.
   * @param r Observation point.
   * @param rp Source point.
   * @return Complex scalar \f$ G_s \f$ (used for near/coincident geometry).
   */
  dcmplx Smoothed(rvec r, rvec rp);

  /**
   * @brief Smoothed scalar kernel \f$ G_s \f$ for a specific translation.
   * @param r Observation point.
   * @param rp Source point.
   * @param t Lattice translation applied to the source.
   * @return Complex scalar \f$ G_s \f$.
   */
  dcmplx SmoothedTranslate(rvec r, rvec rp, rvec t);

  /**
   * @brief Source-gradient \f$\nabla^{\prime} G(\mathbf r,\mathbf r^{\prime})\f$.
   * @param r Observation point.
   * @param rp Source point.
   * @return Complex 3-vector gradient.
   */
  cvec Gradient(rvec r, rvec rp);

  /**
   * @brief Source-gradient for a specific translation.
   * @param r Observation point.
   * @param rp Source point.
   * @param t Lattice translation applied to the source.
   * @return Complex 3-vector gradient.
   */
  cvec GradientTranslate(rvec r, rvec rp, rvec t);

  /**
   * @brief Smoothed source-gradient \f$\nabla^{\prime} G_s\f$.
   * @param r Observation point.
   * @param rp Source point.
   * @return Complex 3-vector gradient (smoothed).
   */
  cvec GradientSmoothed(rvec r, rvec rp);

  /**
   * @brief Smoothed source-gradient for a specific translation.
   * @param r Observation point.
   * @param rp Source point.
   * @param t Lattice translation applied to the source.
   * @return Complex 3-vector gradient (smoothed).
   */
  cvec GradientSmoothedTranslate(rvec r, rvec rp, rvec t);

  /**
   * @brief Evaluate \f$\{G,\ \nabla^{\prime} G\}\f$ for a translation in one call.
   * @param r Observation point.
   * @param rp Source point.
   * @param t Lattice translation applied to the source.
   * @return Pair \f$(G,\ \nabla^{\prime} G)\f$.
   */
  std::pair<dcmplx, cvec> EvaluateandGradientTranslate(rvec r, rvec rp, rvec t);

  /**
   * @brief Evaluate \f$\{G_s,\ \nabla^{\prime} G_s\}\f$ (smoothed) for a translation in one call.
   * @param r Observation point.
   * @param rp Source point.
   * @param t Lattice translation applied to the source.
   * @return Pair \f$(G_s,\ \nabla^{\prime} G_s)\f$.
   */
  std::pair<dcmplx, cvec> EvaluateandGradientSmoothedTranslate(rvec r, rvec rp, rvec t);

  // ---------------------------------------------------------------------------
  // RWG‑weighted integrals
  // ---------------------------------------------------------------------------
  /**
   * @brief D-operator matrix element 
   * \f$\langle \mathbf{f}(\mathbf{r}),\,\overline{\overline G},\,\mathbf{f}(\mathbf{r}^{\prime})\rangle\f$.
   * @param f Testing RWG function (on observation triangle).
   * @param fp Basis RWG function (on source triangle).
   * @return Complex scalar entry of \f$\mathbf D\f$.
   * @note Not used for periodic homogeneous case; exits.
   */
  dcmplx IntegrateD(RWGFun* f, RWGFun* fp);

  /**
   * @brief K-operator matrix element 
   * \f$\langle \mathbf{f}(\mathbf{r}),\,\nabla \times \overline{\overline G},\,\mathbf{f}(\mathbf{r}^{\prime})\rangle\f$.
   * @param f Testing RWG function (on observation triangle).
   * @param fp Basis RWG function (on source triangle).
   * @return Complex scalar entry of \f$\mathbf K\f$.
   * @note Not used for periodic homogeneous case; exits.
   */
  dcmplx IntegrateK(RWGFun* f, RWGFun* fp);

  /**
   * @brief @f$\boldsymbol \delta@f$ term at an observation point.
   * @param r Observation point \f$\mathbf r\f$.
   * @param fp Source RWG function.
   * @return Complex 3-vector contribution \f$\boldsymbol \delta(\mathbf r)\f$.
   */
  cvec IntegrateDelta(rvec r, RWGFun* fp);

  /**
   * @brief @f$\boldsymbol \kappa@f$ term at an observation point.
   * @param r Observation point \f$\mathbf r\f$.
   * @param fp Source RWG function.
   * @return Complex 3-vector contribution \f$\boldsymbol \kappa(\mathbf r)\f$.
   */
  cvec IntegrateKappa(rvec r, RWGFun* fp);

  /**
   * @brief Pair of D elements with swapped triangle roles.
   * @param f Testing RWG function.
   * @param fp Source RWG function.
   * @return \f$( D(f,fp),\, D(fp,f) )\f$.
   */
  std::pair<dcmplx, dcmplx> PairD(RWGFun* f, RWGFun* fp);

  /**
   * @brief Pair of K elements with swapped triangle roles.
   * @param f Testing RWG function.
   * @param fp Source RWG function.
   * @return \f$( K(f,fp),\, K(fp,f) )\f$.
   */
  std::pair<dcmplx, dcmplx> PairK(RWGFun* f, RWGFun* fp);

  // ---------------------------------------------------------------------------
  // Batched helpers
  // ---------------------------------------------------------------------------
  /**
   * @brief Compute \f$\mathbf D\f$ and \f$\mathbf K\f$ for all RWG pairs on two triangles.
   * @param fvec RWGs living on the observation triangle.
   * @param fpvec RWGs living on the source triangle.
   * @param DK Output matrix of pairs \f$(D,K)\f$ indexed by \c [i][j].
   */
  virtual void SameTriDK(const std::vector<RWGFun*>& fvec, 
                         const std::vector<RWGFun*>& fpvec,
                         std::vector<std::vector<std::pair<dcmplx, dcmplx> > >& DK);
  
  /**
   * @brief Compute \f$(\boldsymbol \delta,\boldsymbol \kappa)\f$ for all RWGs on a source triangle.
   * @param r Observation point.
   * @param fpvec RWGs living on the source triangle.
   * @param DeltaKappa Output vector of pairs \f$(\boldsymbol \delta,\boldsymbol \kappa)\f$.
   */
  virtual void SameTriDeltaKappa(rvec r, const std::vector<RWGFun*>& fpvec,
                                 std::vector<std::pair<cvec, cvec> >& DeltaKappa);
};

#endif
