/**
 * @file GreenFHom3DPer.cpp
 * @brief Implementation of the periodic homogeneous 3D Green's function base.
 *
 * Provides only the constructor wiring, a basic Assign() helper used by
 * near/far tests, and a setter for the number of evaluation terms. Concrete
 * periodic kernels derive from this class and implement the actual summations
 * and the PrimitiveCell() mapping.
 */

#include "GreenFHom3DPer.h"

// -----------------------------
// Construction & lifetime
// -----------------------------
GreenFHom3DPer::GreenFHom3DPer(dcmplx wvl, dcmplx eps, rvec kin)
    : GreenFHom3D(wvl, eps), k(kin), evaluationTerms(0) {}

GreenFHom3DPer::GreenFHom3DPer(dcmplx wvl, dcmplx eps, dcmplx mu, rvec kin)
    : GreenFHom3D(wvl, eps, mu), k(kin), evaluationTerms(0) {}

GreenFHom3DPer::~GreenFHom3DPer() {}

// -----------------------------
// Thresholding helper
// -----------------------------
void GreenFHom3DPer::Assign(rvec r, rvec rp, rvec l,
                            std::pair<bool, rvec>& out) {
  if (!GreenFHom3D::AboveThreshold(r - rp, l)) {
    out.first = false;
    out.second = l;
  }
}

// -----------------------------
// Controls
// -----------------------------
void GreenFHom3DPer::EvaluationTerms(int n) {
  std::cout << "Number of evaluation terms: " << n << std::endl;
  evaluationTerms = n;
}
