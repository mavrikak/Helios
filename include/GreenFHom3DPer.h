/**
 * @file GreenFHom3DPer.h
 * @brief Periodic homogeneous 3D Green's function (abstract base).
 *
 * This interface extends #GreenFHom3D to handle **Bloch‑periodic** media in a
 * homogeneous background. Concrete subclasses implement the lattice geometry
 * (primitive cell mapping) and the periodic summation strategy (e.g. Ewald
 * splitting in real/reciprocal space).
 *
 * Key ideas:
 * - Fields satisfy Floquet–Bloch condition with wavevector **k** (a.k.a. Bloch
 * vector). Periodic images are accumulated with phase factors exp(i k·R).
 * - Convergence is typically accelerated with **Ewald** splitting;
 * `evaluationTerms` and `e` control how many terms and the splitting strength.
 * - `PrimitiveCell()` maps any lattice vector to a canonical representative of
 * the primitive cell.
 *
 * @warning This class is **abstract**: `PrimitiveCell()` is pure virtual.
 */

#ifndef GREENFHOM3DPER_H
#define GREENFHOM3DPER_H

#include <string>
#include "GreenFHom3D.h"
#include "globals.h"
#include "mathFunctions.h"

/// \ingroup greenFunction
/**
 * @class GreenFHom3DPer
 * @brief Base class for homogeneous 3D Green's functions with lattice periodicity.
 * 
 * @details Stores common periodic‑summation state (Bloch vector **k**, unit‑cell volume,
 * Ewald/summation controls) and provides a small helper (`Assign`) used in
 * near/far thresholding when deciding which terms to evaluate.
 */
class GreenFHom3DPer : public GreenFHom3D {
 protected:
  // ---------------------------------------------------------------------------
  // Periodic state
  // ---------------------------------------------------------------------------
  rvec k;   //!< Bloch vector.
  double V; //!< Unit‑cell volume.

  // ---------------------------------------------------------------------------
  // Convergence / thresholding controls
  // ---------------------------------------------------------------------------
  std::pair<bool, rvec> aboveThreshold; //!< Cached result of a near/far test.
  int evaluationTerms;                  //!< Number of periodic terms to evaluate.
  double e;                             //!< Ewald splitting parameter.

  /**
   * @brief Map a lattice vector to the primitive cell (subclass must implement).
   * @param target Any vector expressed in the lattice basis.
   * @return The representative inside the primitive cell.
   */
  virtual rvec PrimitiveCell(rvec& target) = 0;
  
  /**
   * @brief Helper used by near/far tests to assign the nearest image vector.
   * @param r Observation point.
   * @param rp Source point.
   * @param testvector Candidate lattice translation vector.
   * @param result Output pair: {isFarEnough, correctedVector}.
   * @details Sets `result.first=false` (i.e., not far enough) and stores the
   * vector to be used when #GreenFHom3D::AboveThreshold fails.
   */
  void Assign(rvec r, rvec rp, rvec testvector, std::pair<bool, rvec>& result);

 public:
  // ---------------------------------------------------------------------------
  // Construction & lifetime
  // ---------------------------------------------------------------------------
  /**
   * @brief Homogeneous, non‑magnetic constructor with Bloch vector.
   * @param wavelength Free‑space wavelength @f$\lambda_0@f$.
   * @param epsilonMedium Relative permittivity @f$\varepsilon_r@f$.
   * @param wavevector Bloch vector @f$\mathbf k@f$.
   */
  GreenFHom3DPer(dcmplx wavelength, dcmplx epsilonMedium, rvec wavevector);

  /**
   * @brief Homogeneous, magnetic constructor with Bloch vector.
   * @param wavelength Free‑space wavelength @f$\lambda_0@f$.
   * @param epsilonMedium Relative permittivity @f$\varepsilon_r@f$.
   * @param muMedium Relative permeability @f$\mu_r@f$.
   * @param wavevector Bloch vector @f$\mathbf k@f$.
   */
  GreenFHom3DPer(dcmplx wavelength, dcmplx epsilonMedium, dcmplx muMedium,
                 rvec wavevector);

  /** @brief Destructor. */
  virtual ~GreenFHom3DPer();
  
  // ---------------------------------------------------------------------------
  // Controls
  // ---------------------------------------------------------------------------
  /**
  * @brief Set the number of periodic evaluation terms to accumulate.
  * @param n Number of periodic evaluation terms.
  * @details Useful for accuracy/performance trade‑offs; echoed to stdout for transparency.
  */
  void EvaluationTerms(int n);
};

#endif
