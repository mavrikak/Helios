/**
 * @file GreenF.h
 * @brief Abstract base class for electromagnetic Green's functions.
 *
 * This class defines the common interface and shared data used by concrete
 * Green's function providers in the codebase (e.g., free‑space, periodic,
 * and layered‑media variants). It exposes pointwise evaluations (scalar,
 * gradient, dyadic), as well as the RWG‑weighted integrals required to
 * assemble MoM/SIE matrices (D/K) and post‑processing vectors (delta/kappa).
 */

#ifndef GREENF_H
#define GREENF_H

#include <blitz/array.h>
#include <cmath>
#include <utility>
#include "GaussQuad.h"
#include "RWGFun.h"
#include "SingSub.h"
#include "globals.h"
#include "mathFunctions.h"

/// \ingroup greenFunction
/**
 * @class GreenF
 * @brief Abstract base for Green's function kernels.
 *
 * Responsibilities:
 * - Store wavelength and wave numbers for homogeneous (k_B) or layered (k_L) media.
 * - Provide a consistent polymorphic API for evaluating kernels and their RWG integrals.
 * - Offer optional acceleration via an erfc lookup table (EvErfc).
 * - Expose accuracy and periodic‑term switches used by some derived classes.
 */
class GreenF {
 protected:
  // ---------------------------------------------------------------------------
  // Medium parameters
  // ---------------------------------------------------------------------------
  dcmplx wavelength;        //!< Wavelength in medium.
  dcmplx k_B;               //!< Homogeneous wavenumber (if applicable).
  std::vector<dcmplx> k_L;  //!< Per‑layer wavenumbers for layered problems.

  // ---------------------------------------------------------------------------
  // erfc() lookup table (optional acceleration)
  // ---------------------------------------------------------------------------
  LookupTableBin* t;          //!< Pointer to error‑function lookup table (owned elsewhere).
  bool thereIsNoLookupTable;  //!< If true, fall back to direct cerfc() evaluation.

  // ---------------------------------------------------------------------------
  // Global switches
  // ---------------------------------------------------------------------------
  static bool NeedsAccurate;  //!< Request high‑accuracy integration where available.
  static int Etm;             //!< Number of periodic evaluation terms (for periodic GFs).

 public:
  // ---------------------------------------------------------------------------
  // Construction & lifetime
  // ---------------------------------------------------------------------------
  /**
   * @brief Constructor for a homogeneous, non-magnetic medium.
   * @param wavelength Wavelength in the medium.
   * @param epsilonMedium Relative permittivity \f$\varepsilon_r\f$.
   */
  GreenF(dcmplx wavelength, dcmplx epsilonMedium);
  
  /**
   * @brief Constructor for a layered, non-magnetic medium.
   * @param wavelength Wavelength in the medium.
   * @param epsilonMedium Vector of layer permittivities \f$\varepsilon_{r,\ell}\f$.
   */
  GreenF(dcmplx wavelength, const std::vector<dcmplx>& epsilonMedium);
  
  /**
   * @brief Constructor for a homogeneous, magnetic medium.
   * @param wavelength Wavelength in the medium.
   * @param epsilonMedium Relative permittivity \f$\varepsilon_r\f$.
   * @param muMedium Relative permeability \f$\mu_r\f$.
   */
  GreenF(dcmplx wavelength, dcmplx epsilonMedium, dcmplx muMedium);
  
  /**
   * @brief Constructor for a layered, magnetic medium.
   * @param wavelength Wavelength in the medium.
   * @param epsilonMedium Vector of layer permittivities \f$\varepsilon_{r,\ell}\f$.
   * @param muMedium Vector of layer permeabilities \f$\mu_{r,\ell}\f$.
   */
  GreenF(dcmplx wavelength, 
         const std::vector<dcmplx>& epsilonMedium, 
         const std::vector<dcmplx>& muMedium);

  /** @brief Destructor. */
  virtual ~GreenF();

  // ---------------------------------------------------------------------------
  // Pointwise kernel evaluation
  // ---------------------------------------------------------------------------
  /**
   * @brief Scalar Green's function \f$ G(\mathbf r,\mathbf r^{\prime}) \f$.
   * @param r Observation point \f$\mathbf r\f$.
   * @param rp Source point \f$\mathbf r^{\prime}\f$.
   * @return Complex scalar \f$ G(\mathbf r,\mathbf r^{\prime})\f$.
   */
  virtual dcmplx Evaluate(rvec r, rvec rp) = 0;

  /**
   * @brief Gradient w.r.t. the source point \f$\nabla^{\prime}G(\mathbf r,\mathbf r^{\prime})\f$.
   * @param r Observation point \f$\mathbf r\f$.
   * @param rp Source point \f$\mathbf r^{\prime}\f$.
   * @return Complex 3-vector \f$\nabla^{\prime} G(\mathbf r,\mathbf r^{\prime})\f$.
   */
  virtual cvec Gradient(rvec r, rvec rp) = 0;

  /**
   * @brief Dyadic Green's function \f$\overline{\overline G}(\mathbf r,\mathbf r^{\prime})\f$.
   * @param r Observation point \f$\mathbf r\f$.
   * @param rp Source point \f$\mathbf r^{\prime}\f$.
   * @return Complex \f$ 3 \times 3 \f$ dyad.
   */
  virtual cdyad EvaluateDyadic(rvec r, rvec rp) = 0;

  // ---------------------------------------------------------------------------
  // RWG‑weighted integrals for MoM/SIE assembly
  // ---------------------------------------------------------------------------
  /**
   * @brief D-operator matrix element 
   * \f$\langle \mathbf{f}(\mathbf{r}),\,\overline{\overline G},\,\mathbf{f}(\mathbf{r}^{\prime})\rangle\f$.
   * @param f Testing RWG function (on observation triangle).
   * @param fp Basis RWG function (on source triangle).
   * @return Complex scalar entry of \f$\mathbf D\f$.
   */
  virtual dcmplx IntegrateD(RWGFun* f, RWGFun* fp) = 0;

  /**
   * @brief K-operator matrix element 
   * \f$\langle \mathbf{f}(\mathbf{r}),\,\nabla \times \overline{\overline G},\,\mathbf{f}(\mathbf{r}^{\prime})\rangle\f$.
   * @param f Testing RWG function (on observation triangle).
   * @param fp Basis RWG function (on source triangle).
   * @return Complex scalar entry of \f$\mathbf K\f$.
   */
  virtual dcmplx IntegrateK(RWGFun* f, RWGFun* fp) = 0;

  // ---------------------------------------------------------------------------
  // Post‑processing (fields from solved currents)
  // ---------------------------------------------------------------------------
  /**
   * @brief @f$\boldsymbol \delta@f$ term at an observation point.
   * @param r Observation point \f$\mathbf r\f$.
   * @param fp Source RWG function.
   * @return Complex 3-vector contribution \f$\boldsymbol \delta(\mathbf r)\f$.
   */
  virtual cvec IntegrateDelta(rvec r, RWGFun* fp) = 0;

  /**
   * @brief @f$\boldsymbol \kappa@f$ term at an observation point.
   * @param r Observation point \f$\mathbf r\f$.
   * @param fp Source RWG function.
   * @return Complex 3-vector contribution \f$\boldsymbol \kappa(\mathbf r)\f$.
   */
  virtual cvec IntegrateKappa(rvec r, RWGFun* fp) = 0;

  /**
   * @brief Pair of D elements with swapped triangle roles.
   * @param f Testing RWG function.
   * @param fp Basis RWG function.
   * @return \f$( D(f,fp),\, D(fp,f) )\f$.
   */
  virtual std::pair<dcmplx, dcmplx> PairD(RWGFun* f, RWGFun* fp) = 0;

  /**
   * @brief Pair of K elements with swapped triangle roles.
   * @param f Testing RWG function.
   * @param fp Basis RWG function.
   * @return \f$( K(f,fp),\, K(fp,f) )\f$.
   */
  virtual std::pair<dcmplx, dcmplx> PairK(RWGFun* f, RWGFun* fp) = 0;

  // ---------------------------------------------------------------------------
  // D and K matrix elements computation
  // ---------------------------------------------------------------------------
  /**
   * @brief Compute \f$\mathbf D\f$ and \f$\mathbf K\f$ for all RWG pairs on two triangles.
   * @param fvec RWGs living on the observation triangle.
   * @param fpvec RWGs living on the source triangle.
   * @param DK Output matrix of pairs \f$(D,K)\f$ indexed by \c [i][j].
   */
  virtual void SameTriDK(const std::vector<RWGFun*>& fvec, 
                         const std::vector<RWGFun*>& fpvec,
                         std::vector<std::vector<std::pair<dcmplx, dcmplx>>>& DK) = 0;
  
  /**
   * @brief Same as SameTriDK plus gradients (if provided by the GF).
   * @param fvec RWGs on the observation triangle.
   * @param fpvec RWGs on the source triangle.
   * @param DK Output matrix of \f$(D,K)\f$ pairs.
   * @param DKG Output matrix of gradient pairs per entry.
   */
  virtual void SameTriDKG(const std::vector<RWGFun*>& fvec, 
                          const std::vector<RWGFun*>& fpvec,
                          std::vector<std::vector<std::pair<dcmplx, dcmplx>>>& DK,
                          std::vector<std::vector<std::pair<cvec, cvec>>>& DKG) = 0;

  /**
   * @brief Compute \f$(\boldsymbol \delta,\boldsymbol \kappa)\f$ for all RWGs on a source triangle.
   * @param r Observation point.
   * @param fpvec RWGs living on the source triangle.
   * @param DeltaKappa Output vector of pairs \f$(\boldsymbol \delta,\boldsymbol \kappa)\f$.
   */
  virtual void SameTriDeltaKappa(rvec r, const std::vector<RWGFun*>& fpvec,
                                 std::vector<std::pair<cvec, cvec>>& DeltaKappa) = 0;

  // ---------------------------------------------------------------------------
  // Utilities
  // ---------------------------------------------------------------------------
  /** @brief Return homogeneous wavevector \f$ k_B \f$ (if defined).
   *  @return Complex wavenumber \f$ k_B \f$. 
   */
  dcmplx WaveVec();

  /**
   * @brief Evaluate \c erfc via the 2D lookup table if available.
   * @param z Complex argument.
   * @return \c erfc(z) value (from table or direct evaluation).
   */
  dcmplx EvErfc(dcmplx z);

  /**
   * @brief Initialize the \c erfc lookup table pointer.
   * @param t Pointer to a binary lookup table (owned elsewhere).
   * @return 0 on success; non-zero otherwise.
   */
  int InitLookupTable(LookupTableBin* t);

  /** @brief Request high‑accuracy evaluation paths (global). */
  static void EnableAccurate();

  /** @brief Query whether high-accuracy evaluation is requested.
   *  @return true if accurate paths are enabled; false otherwise. 
   */
  static bool RequiresAccurate();

  /**
   * @brief Set number of periodic evaluation terms (for periodic GFs).
   * @param n Non-negative number of terms to include.
   */
  static void AssignEtm(int n);

  /**
   * @brief Get layer wavenumber by index.
   * @param index Layer index (0-based).
   * @return Complex wavenumber \f$ k_{L,\mathrm{index}} \f$.
   * @throws std::out_of_range if @p index is invalid.
   */
  dcmplx GetLayerWavenumber(size_t index) const;

  /**
   * @brief Get all layer wavenumbers.
   * @return Vector of complex wavenumbers \f$\{k_{L,\ell}\}\f$.
   */
  std::vector<dcmplx> GetLayersWavenumbers() const;
};

/**
 * @ingroup greenFunction
 * @brief Complex square root with a branch cut along the negative real axis.
 * @param arg Argument of square root.
 * @return Result of square root operation.
 * @details Ensures continuity for purely real negative inputs by nudging the
 * imaginary part to zero before calling std::sqrt.
 */
dcmplx csqrt(dcmplx arg);

#endif
