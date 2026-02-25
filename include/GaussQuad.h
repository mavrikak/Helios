/**
 * @file GaussQuad.h
 * @brief Gaussian quadrature helpers for Green's functions and incident fields.
 *
 * The templated class #GaussQuad<Type> performs first‑order Gaussian
 * quadrature on triangular surface elements for two broad families of
 * integrands:
 *
 * 1) **Green's function kernels** used in SIE matrix assembly.
 * 2) **Incident fields** used to build right‑hand‑side vectors.
 *
 * The interface is intentionally minimal and relies on the following contracts
 * with the template parameter @p Type:
 * - For **Green's integrals**: @p Type must expose dyadic or scalar Green's
 * evaluation methods with signature
 * - `dcmplx Type::(rvec r, rvec rp)` for scalar kernels, and/or
 * - `cvec Type::(rvec r, rvec rp)` for vector kernels,
 * and a `WaveVec()` method returning @f$ k @f$.
 * - For **incident integrals**: @p Type must expose
 * - `cvec Type::(rvec r)` (single field), or
 * - `std::vector<cvec> Type::(rvec r)` returning `{E, H}` at @p r.
 *
 * Integration accuracy is controlled by Dunavant rules (see `quadrature.h`).
 */

#ifndef GAUSSQUAD_H
#define GAUSSQUAD_H

#include <cmath>
#include <iostream>
#include "RWGFun.h"
#include "globals.h"
#include "quadrature.h"

/// \ingroup quadrature
/**
 * @class GaussQuad
 * @brief Gaussian-quadrature helpers for Green's functions and incident fields.
 * @tparam Type either a Green's function provider or an incident‑field class
 * exposing the member functions described in the file header.
 *
 * @note The implementation uses **centroid evaluation** (first order) for the
 * inner integrals `I1`, `I2`, `I4`, `I5`, which is consistent with the
 * low‑order RWG discretization and the JOSA B Vol 26 No 4 2009 first‑order formulas
 * referenced in the code comments below.
 */
template <class Type>
class GaussQuad {
 private:
  // ---------------------------------------------------------------------------
  // Inner integrals (centroid evaluation) used by the first‑order formulas
  // ---------------------------------------------------------------------------
  /** @brief I1: integral of scalar kernel times RWG prefactor over the test triangle. */
  dcmplx I1(dcmplx grrp, RWGFun* fp);
  /** @brief I2: vector form with the RWG basis evaluated at the centroid. */
  cvec I2(dcmplx grrp, RWGFun* fp);
  /** @brief I4: vector kernel crossed with RWG at the centroid. */
  cvec I4(cvec grrp, RWGFun* fp);

 public:
  // ---------------------------------------------------------------------------
  // Green's integrals on one or two triangles
  // ---------------------------------------------------------------------------
  /**
   * @brief I1 with on‑the‑fly kernel evaluation at (r,rp) using @p gF.
   * @param gF Green's function object.
   * @param EvaluationFunction Pointer to scalar kernel member of @p gF.
   * @param r Observation point (usually centroid of the test triangle).
   * @param fp Source RWG living on the source triangle.
   * @return Complex scalar integral value.
   */
  dcmplx I1(Type* gF, dcmplx (Type::*)(rvec r, rvec rp), rvec r, RWGFun* fp);

  /**
   * @brief I2 with on‑the‑fly scalar kernel evaluation.
   * @param gF Green's function object.
   * @param EvaluationFunction Pointer to scalar kernel member of @p gF.
   * @param r Observation point (usually centroid of the test triangle).
   * @param fp Source RWG living on the source triangle.
   * @return Complex vector integral value.
   */
  cvec I2(Type* gF, dcmplx (Type::*)(rvec r, rvec rp), rvec r, RWGFun* fp);

  /**
   * @brief I4 with on‑the‑fly vector kernel evaluation.
   * @param gF Green's function object.
   * @param EvaluationFunction Pointer to scalar kernel member of @p gF.
   * @param r Observation point (usually centroid of the test triangle).
   * @param fp Source RWG living on the source triangle.
   * @return Complex vector integral value.
   */
  cvec I4(Type* gF, cvec (Type::*)(rvec r, rvec rp), rvec r, RWGFun* fp);

  /**
   * @brief I5 variant: vector kernel weighted by RWG prefactor (centroid).
   * @param gF Green's function object.
   * @param EvaluationFunction Pointer to scalar kernel member of @p gF.
   * @param r Observation point (usually centroid of the test triangle).
   * @param fp Source RWG living on the source triangle.
   * @return Complex vector integral value.
   */
  cvec I5(Type* gF, cvec (Type::*)(rvec r, rvec rp), rvec r, RWGFun* fp);
  
  /**
   * @brief Two‑triangle integral for the **D** operator (Eq. (27), JOSA B 26(4), 2009).
   * @param gF Green's function object.
   * @param EvaluationFunction Pointer to scalar kernel member of @p gF.
   * @param f Observation RWG living on the observation triangle.
   * @param fp Source RWG living on the source triangle.
   * @return Complex scalar value.
   * @details Uses first‑order centroid approximations of inner integrals and
   * the wavenumber @f$ k @f$ provided by `Type::WaveVec()`.
   */
  dcmplx IntegrateD(Type* gF, dcmplx (Type::*)(rvec r, rvec rp), RWGFun* f,
                    RWGFun* fp);
  
  /**
   * @brief Two‑triangle integral for the **K** operator (Eq. (29), JOSA B 26(4), 2009).
   * @param gF Green's function object.
   * @param EvaluationFunction Pointer to scalar kernel member of @p gF.
   * @param f Observation RWG living on the observation triangle.
   * @param fp Source RWG living on the source triangle.
   * @return Complex scalar value.
   * @details Uses first‑order centroid approximations of inner integrals and
   * the wavenumber @f$ k @f$ provided by `Type::WaveVec()`.
   */
  dcmplx IntegrateK(Type* gF, cvec (Type::*)(rvec r, rvec rp), RWGFun* f,
                    RWGFun* fp);

  // ---------------------------------------------------------------------------
  // Incident‑field integrals (RHS)
  // ---------------------------------------------------------------------------
  /**
   * @brief Integrate on triangle E_inc using centroid evaluation.
   * @param gF Incident field provider.
   * @param EvaluationFunction Pointer to `cvec Type::(rvec)`.
   * @param f Test RWG function.
   * @return Complex scalar value.
   */
  dcmplx Integrate(Type* gF, cvec (Type::*)(rvec r), RWGFun* f);

  /**
   * @brief Integrate on triangle @f$ \{ \mathbf{E}^{inc}, \mathbf{H}^{inc} \} @f$ with Dunavant quadrature.
   * @param gF Incident field provider.
   * @param EvaluationFunction Pointer to `std::vector<cvec> Type::(rvec)` that returns `{E,H}`.
   * @param f Test RWG function.
   * @param NeedsAccurate When true, use a higher‑order Dunavant rule (N5); otherwise N1.
   * @return Vector @f$ \{\iint \mathbf{f} \cdot \mathbf{E}^{inc} dS, \iint \mathbf{f} \cdot \mathbf{H}^{inc} dS \} @f$ in this order.
   */
  std::vector<dcmplx> Integrate(
        Type* gF, std::vector<cvec> (Type::*EvaluationFunction)(rvec r), 
        RWGFun* f, bool NeedsAccurate);
};

// ===================== Implementation (inline) =====================

template <class Type>
dcmplx GaussQuad<Type>::I1(Type* gF,
                                dcmplx (Type::*EvaluationFunction)(rvec, rvec),
                                rvec r, RWGFun* fp) {
  rvec rp(fp->TrianglePtr()->Center());
  dcmplx grrp((*gF.*EvaluationFunction)(r, rp));
  double Areap(fp->TrianglePtr()->Area());
  return Areap * fp->RWGPreFactor() * grrp;
}

template <class Type>
cvec GaussQuad<Type>::I2(Type* gF,
                              dcmplx (Type::*EvaluationFunction)(rvec,
                                                                      rvec),
                              rvec r, RWGFun* fp) {
  rvec rp(fp->TrianglePtr()->Center());
  dcmplx grrp((*gF.*EvaluationFunction)(r, rp));
  double Areap(fp->TrianglePtr()->Area());
  return Areap * fp->Evaluate(rp) * grrp;
}

template <class Type>
cvec GaussQuad<Type>::I4(Type* gF,
                              cvec (Type::*EvaluationFunction)(rvec, rvec),
                              rvec r, RWGFun* fp) {
  rvec rp(fp->TrianglePtr()->Center());
  cvec grrp((*gF.*EvaluationFunction)(r, rp));
  cvec fprp(fp->Evaluate(rp));
  double Areap(fp->TrianglePtr()->Area());
  return Areap * cross(grrp, fprp);
}

template <class Type>
cvec GaussQuad<Type>::I5(Type* gF,
                              cvec (Type::*EvaluationFunction)(rvec, rvec),
                              rvec r, RWGFun* fp) {
  rvec rp(fp->TrianglePtr()->Center());
  cvec grrp((*gF.*EvaluationFunction)(r, rp));
  double Areap(fp->TrianglePtr()->Area());
  return Areap * fp->RWGPreFactor() * grrp;
}

// Returns evaluation of equation (27) in JOSA B Vol 26 No 4 2009 at first order.
template <class Type>
dcmplx GaussQuad<Type>::IntegrateD(
    Type* gF, dcmplx (Type::*EvaluationFunction)(rvec, rvec),
    RWGFun* f, RWGFun* fp) {
  rvec r(f->TrianglePtr()->Center());
  rvec rp(fp->TrianglePtr()->Center());
  double Area(f->TrianglePtr()->Area());
  dcmplx grrp((*gF.*EvaluationFunction)(r, rp));
  dcmplx ki(gF->WaveVec());
  return Area * (dot(f->Evaluate(r), I2(grrp, fp)) -
                 1. / (ki * ki) * f->RWGPreFactor() * I1(grrp, fp));
}

// Returns evaluation of equation (29) in JOSA B Vol 26 No 4 2009 at first order.
template <class Type>
dcmplx GaussQuad<Type>::IntegrateK(
    Type* gF, cvec (Type::*EvaluationFunction)(rvec r, rvec rp),
    RWGFun* f, RWGFun* fp) {
  rvec r(f->TrianglePtr()->Center());
  rvec rp(fp->TrianglePtr()->Center());
  double Area(f->TrianglePtr()->Area());
  cvec grrp((*gF.*EvaluationFunction)(r, rp));
  cvec fr(f->Evaluate(r));
  return Area * dot(fr, I4(grrp, fp));
}

template <class Type>
dcmplx GaussQuad<Type>::Integrate(
    Type* inF, cvec (Type::*EvaluationFunction)(rvec r),
    RWGFun* f) {
  double Area((*f).TrianglePtr()->Area());
  rvec r(0., 0., 0.);
  for (int i(0); i <= 2; ++i) {
    r += f->TrianglePtr()->Node(i);
  }
  r /= 3.;
  cvec rf(f->Evaluate(r));
  dcmplx out(Area * dot(rf, (*inF.*EvaluationFunction)(r)));
  return out;
}

template <class Type>
std::vector<dcmplx> GaussQuad<Type>::Integrate(
    Type* inF, std::vector<cvec> (Type::*EvaluationFunction)(rvec r),
    RWGFun* f, bool NeedsAccurate) {
  double Area((*f).TrianglePtr()->Area());
  rvec r[3];
  for (int i = 0; i < 3; i++) {
    r[i] = f->TrianglePtr()->Node(i);
  }
  int Ndunavant;
  double(*xdunavant)[3];
  double* wdunavant;
  Ndunavant = NeedsAccurate ? dunavant::N5 : dunavant::N1;
  xdunavant = NeedsAccurate ? dunavant::x5 : dunavant::x1;
  wdunavant = NeedsAccurate ? dunavant::w5 : dunavant::w1;
  std::vector<rvec> svec(Ndunavant);
  for (int i = 0; i < Ndunavant; ++i) {
    svec[i] = r[0] * xdunavant[i][0] + r[1] * xdunavant[i][1] +
              r[2] * xdunavant[i][2];
  }
  dcmplx outE(0.), outH(0.);
  for (int i = 0; i < Ndunavant; ++i){
    double dA = Area * wdunavant[i];
    cvec rf(f->Evaluate(svec[i]));
    outE += dA * dot(rf, ((*inF.*EvaluationFunction)(svec[i])).at(0));
    outH += dA * dot(rf, ((*inF.*EvaluationFunction)(svec[i])).at(1));
  }
  return {outE, outH};
}

template <class Type>
dcmplx GaussQuad<Type>::I1(dcmplx grrp, RWGFun* fp) {
  rvec rp(fp->TrianglePtr()->Center());
  double Areap(fp->TrianglePtr()->Area());
  return Areap * fp->RWGPreFactor() * grrp;
}

template <class Type>
cvec GaussQuad<Type>::I2(dcmplx grrp, RWGFun* fp) {
  rvec rp(fp->TrianglePtr()->Center());
  double Areap(fp->TrianglePtr()->Area());
  return Areap * fp->Evaluate(rp) * grrp;
}

template <class Type>
cvec GaussQuad<Type>::I4(cvec grrp, RWGFun* fp) {
  rvec rp(fp->TrianglePtr()->Center());
  cvec fprp(fp->Evaluate(rp));
  double Areap(fp->TrianglePtr()->Area());
  return Areap * cross(grrp, fprp);
}
#endif
