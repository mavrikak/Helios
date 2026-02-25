/**
 * @file GreenFHom3DPer2D.cpp
 * @brief Implementation of the 2D-periodic homogeneous 3D Green's function.
 *
 * This file provides the Bloch-periodic (2D lattice) free-space Green's
 * function and its utilities:
 * - Ewald-split evaluation of the scalar kernel G and source-gradient grad'G,
 * including smoothed variants used in singularity subtraction.
 * - RWG-weighted point–triangle integrals for vectors used in the SIE RHS 
 * and post-processing.
 * - Batched same-triangle helpers for efficient assembly.
 *
 * Important: This is a specialization of GreenFHom3DPer; triangle–triangle
 * assembly of D/K in the periodic free-space case is disabled in favor of
 * layered kernels. Point–triangle routines are used for RHS/fields.
 */

#include "GreenFHom3DPer2D.h"
#include <fstream>
#include <iostream>
#include "quadrature.h"

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------
GreenFHom3DPer2D::GreenFHom3DPer2D(dcmplx wvl, dcmplx eps, rvec kin, rvec lat1,
                                   rvec lat2)
    : GreenFHom3DPer(wvl, eps, kin), a1(lat1), a2(lat2) {
  Init();
}

GreenFHom3DPer2D::GreenFHom3DPer2D(dcmplx wvl, dcmplx eps, dcmplx mu, rvec kin,
                                   rvec lat1, rvec lat2)
    : GreenFHom3DPer(wvl, eps, mu, kin), a1(lat1), a2(lat2) {
  Init();
}

// ---------------------------------------------------------------------------
// Initialization
// ---------------------------------------------------------------------------
// Common initialization: lattice metrics, reciprocal lattice, k‑projection, Ewald parameter.
void GreenFHom3DPer2D::Init() {
  thereIsNoLookupTable = 1;

  a1a1 = dot(a1, a1);
  a2a2 = dot(a2, a2);
  norma1 = sqrt(a1a1);
  norma2 = sqrt(a2a2);
  a1a2 = dot(a1, a2);
  crossFactor = (a1a1 * a2a2) - a1a2 * a1a2;
  perp = cross(a1, a2);
  perp /= sqrt(dot(perp, perp));

  Reciprocal();
  ToPlane();

  e = PI / sqrt(V);
  EvaluationTerms(Etm);
}

// ---------------------------------------------------------------------------
// Lattice helpers
// ---------------------------------------------------------------------------
// Project Bloch vector onto the lattice plane (discard normal component).
void GreenFHom3DPer2D::ToPlane() {
  double ka1(dot(k, a1));
  double ka2(dot(k, a2));
  double k1((ka1 * a2a2 - ka2 * a1a2) / crossFactor);
  double k2((ka2 * a1a1 - ka1 * a1a2) / crossFactor);
  k = k1 * a1 + k2 * a2;
}

// Build reciprocal lattice vectors (b1, b2) and unit‑cell volume V.
void GreenFHom3DPer2D::Reciprocal() {
  rvec c(cross(a1, a2));
  double norm2(dot(c, c));
  b1 = (2 * PI / norm2) * cross(a2, c);
  b2 = (-2 * PI / norm2) * cross(a1, c);
  V = sqrt(norm2);
}

// Map a displacement to the primitive cell; returns the translation used.
rvec GreenFHom3DPer2D::PrimitiveCell(rvec& target) {
  double ta1(dot(target, a1));
  double ta2(dot(target, a2));
  rvec initial(target);
  double t1((ta1 * a2a2 - ta2 * a1a2) / crossFactor);
  double t2((ta2 * a1a1 - ta1 * a1a2) / crossFactor);

  const double thresh= 0.5+1e-4;
  while ((t1 > thresh || t1 < -thresh)) {
    if (t1 > thresh) {
      target -= a1;
    }
    if (t1 < -thresh) {
      target += a1;
    }
    ta1 = dot(target, a1);
    ta2 = dot(target, a2);
    t1 = (ta1 * a2a2 - ta2 * a1a2) / crossFactor;
  }
  while ((t2 > thresh || t2 < -thresh)) {
    if (t2 > thresh) {
      target -= a2;
    }
    if (t2 < -thresh) {
      target += a2;
    }
    ta1 = dot(target, a1);
    ta2 = dot(target, a2);
    t2 = (ta2 * a1a1 - ta1 * a1a2) / crossFactor;
  }
  return initial - target;
}

// ---------------------------------------------------------------------------
// Thresholds & translations
// ---------------------------------------------------------------------------
//Far/near heuristic for translated triangle pairs (centroid distance vs. size).
bool GreenFHom3DPer2D::AboveThresholdTranslate(Triangle* T, Triangle* Tp,
                                               rvec t) {
  rvec dr(T->Center() - Tp->Center() - t);
  double d = sqrt(dot(dr, dr));
  double psum = T->Perimeter() + Tp->Perimeter();
  return d > psum / 2.;
}

// ---------------------------------------------------------------------------
// Ewald split — spatial (real‑space images)
// ---------------------------------------------------------------------------
// Spatial part of Ewald sum for scalar kernel.
dcmplx GreenFHom3DPer2D::EwaldSpatial(rvec r) {
  dcmplx out(0.);
  double ka1(dot(k, a1));
  double ka2(dot(k, a2));
  for (int i(-evaluationTerms); i <= evaluationTerms; ++i) {
    for (int j(-evaluationTerms); j <= evaluationTerms; ++j) {
      rvec dif(r - i * a1 - j * a2);
      double dis(sqrt(dot(dif, dif)));
      dcmplx sum(exp(-I * k_B * dis) * EvErfc(dis * e - I * k_B / (2 * e)) +
                 exp(I * k_B * dis) * EvErfc(dis * e + I * k_B / (2 * e)));
      out += exp(i * ka1 * I + j * ka2 * I) * sum / dis;
    }
  }
  out /= 8 * PI;
  return out;
}

// Spatial part of Ewald sum for source‑gradient grad'G.
cvec GreenFHom3DPer2D::GradientEwaldSpatial(rvec r) {
  cvec out(0., 0., 0.);
  double ka1(dot(k, a1));
  double ka2(dot(k, a2));
  double spi(sqrt(PI));
  for (int i(-evaluationTerms); i <= evaluationTerms; ++i) {
    for (int j(-evaluationTerms); j <= evaluationTerms; ++j) {
      rvec dif(r - i * a1 - j * a2);
      double dis(sqrt(dot(dif, dif)));
      dcmplx sum(exp(-I * k_B * dis) * (-1. - I * k_B * dis) *
                     EvErfc(dis * e - I * k_B / (2 * e)) +
                 exp(I * k_B * dis) * (-1. + I * k_B * dis) *
                     EvErfc(dis * e + I * k_B / (2 * e)) -
                 dis * 4 * e *
                     exp(-dis * dis * e * e + k_B * k_B / (4 * e * e)) / spi);
      out += exp(i * ka1 * I + j * ka2 * I) * sum * dif / (dis * dis * dis);
    }
  }
  out /= -8 * PI;
  return out;
}

// Spatial part returning {G, grad'G} at once.
std::pair<dcmplx, cvec> GreenFHom3DPer2D::EwaldAndGradientSpatial(rvec r) {
  std::pair<dcmplx, cvec> out(0., cvec(0., 0., 0.));
  double ka1(dot(k, a1));
  double ka2(dot(k, a2));
  double spi(sqrt(PI));
  for (int i(-evaluationTerms); i <= evaluationTerms; ++i) {
    for (int j(-evaluationTerms); j <= evaluationTerms; ++j) {
      rvec dif(r - i * a1 - j * a2);
      double dis(sqrt(dot(dif, dif)));
      dcmplx term1 = exp(-I * k_B * dis) * EvErfc(dis * e - I * k_B / (2 * e));
      dcmplx term2 = exp(I * k_B * dis) * EvErfc(dis * e + I * k_B / (2 * e));
      dcmplx prefactor = exp(i * ka1 * I + j * ka2 * I);
      out.first += prefactor * (term1 + term2) / dis;
      dcmplx sum((-1. - I * k_B * dis) * term1 + (-1. + I * k_B * dis) * term2 -
                 dis * 4 * e *
                     exp(-dis * dis * e * e + k_B * k_B / (4 * e * e)) / spi);
      out.second += prefactor * sum * dif / (dis * dis * dis);
    }
  }
  out.first /= 8 * PI;
  out.second /= -8 * PI;
  return out;
}

// Smoothed spatial part of Ewald sum for scalar kernel (near singularity).
dcmplx GreenFHom3DPer2D::EwaldSpatialSmoothed(rvec r) {
  dcmplx id(1., 0.);
  dcmplx out(0.);
  double ka1(dot(k, a1));
  double ka2(dot(k, a2));
  double eps(1.e-10);
  for (int i(-evaluationTerms); i <= evaluationTerms; ++i) {
    for (int j(-evaluationTerms); j <= evaluationTerms; ++j) {
      rvec dif(r - i * a1 - j * a2);
      double dis(sqrt(dot(dif, dif)));
      if (!((i == 0) && (j == 0))) {  // Regular formula
        dcmplx sum(exp(-I * k_B * dis) * EvErfc(dis * e - I * k_B / (2 * e)) +
                   exp(I * k_B * dis) * EvErfc(dis * e + I * k_B / (2 * e)));
        out += exp(i * ka1 * I + j * ka2 * I) * sum / dis;
      } else {            // On the singularity
        if (dis > eps) {  // Removes singular term numerically
          dcmplx sum(exp(-I * k_B * dis) * EvErfc(dis * e - I * k_B / (2 * e)) +
                     exp(I * k_B * dis) * EvErfc(dis * e + I * k_B / (2 * e)));
          out += (sum - 2. * id) / dis + k_B * k_B * dis;
        } else {  // Takes limit if exactly zero
          out += 2. * I * k_B * (EvErfc(I * k_B / (2 * e)) - id) -
                 4. * e / sqrt(PI) * exp(k_B * k_B / (4. * e * e));
        }
      }
    }
  }
  out /= 8 * PI;
  return out;
}

// Smoothed spatial part of Ewald sum for source‑gradient.
cvec GreenFHom3DPer2D::GradientEwaldSpatialSmoothed(rvec r) {
  dcmplx id(1., 0.);
  cvec out(0., 0., 0.);
  double ka1(dot(k, a1));
  double ka2(dot(k, a2));
  double spi(sqrt(PI));
  double eps(1.e-10);
  for (int i(-evaluationTerms); i <= evaluationTerms; ++i) {
    for (int j(-evaluationTerms); j <= evaluationTerms; ++j) {
      rvec dif(r - i * a1 - j * a2);
      double dis(sqrt(dot(dif, dif)));
      if (!((i == 0) && (j == 0))) {  // Regular formula
        dcmplx sum(exp(-I * k_B * dis) * (-1. - I * k_B * dis) *
                       EvErfc(dis * e - I * k_B / (2 * e)) +
                   exp(I * k_B * dis) * (-1. + I * k_B * dis) *
                       EvErfc(dis * e + I * k_B / (2 * e)) -
                   dis * 4 * e *
                       exp(-dis * dis * e * e + k_B * k_B / (4 * e * e)) / spi);
        out += exp(i * ka1 * I + j * ka2 * I) * sum * dif / (dis * dis * dis);
      } else {            // On the singularity
        if (dis > eps) {  // Removes singular term numerically
          dcmplx sum(exp(-I * k_B * dis) * (-1. - I * k_B * dis) *
                         EvErfc(dis * e - I * k_B / (2 * e)) +
                     exp(I * k_B * dis) * (-1. + I * k_B * dis) *
                         EvErfc(dis * e + I * k_B / (2 * e)) -
                     dis * 4 * e *
                         exp(-dis * dis * e * e + k_B * k_B / (4 * e * e)) /
                         spi);
          out +=
              (sum + 2. * id) * dif / (dis * dis * dis) + k_B * k_B * dif / dis;
        } else {  // Takes limit if exactly zero
          cvec radial(1., 1., 1.);
          out += k_B * k_B * radial / sqrt(3.);
        }
      }
    }
  }
  out /= -8 * PI;
  return out;
}

// Spatial part returning {G, grad'G} with smoothing.
std::pair<dcmplx, cvec> GreenFHom3DPer2D::EwaldAndGradientSpatialSmoothed(
    rvec r) {
  std::pair<dcmplx, cvec> out(0., cvec(0., 0., 0.));
  dcmplx id(1., 0.);
  double ka1(dot(k, a1));
  double ka2(dot(k, a2));
  double spi(sqrt(PI));
  double eps(1.e-10);
  for (int i(-evaluationTerms); i <= evaluationTerms; ++i) {
    for (int j(-evaluationTerms); j <= evaluationTerms; ++j) {
      rvec dif(r - i * a1 - j * a2);
      double dis(sqrt(dot(dif, dif)));
      if (i != 0 || j != 0) {  // Regular formula
        dcmplx term1 =
            exp(-I * k_B * dis) * EvErfc(dis * e - I * k_B / (2 * e));
        dcmplx term2 = exp(I * k_B * dis) * EvErfc(dis * e + I * k_B / (2 * e));
        dcmplx prefactor = exp(i * ka1 * I + j * ka2 * I);
        out.first += prefactor * (term1 + term2) / dis;
        dcmplx sum((-1. - I * k_B * dis) * term1 +
                   (-1. + I * k_B * dis) * term2 -
                   dis * 4 * e *
                       exp(-dis * dis * e * e + k_B * k_B / (4 * e * e)) / spi);
        out.second += prefactor * sum * dif / (dis * dis * dis);
      } else if (dis > eps) {  // Removes singular term numerically
        dcmplx term1 =
            exp(-I * k_B * dis) * EvErfc(dis * e - I * k_B / (2 * e));
        dcmplx term2 = exp(I * k_B * dis) * EvErfc(dis * e + I * k_B / (2 * e));
        out.first += (term1 + term2 - 2. * id) / dis + k_B * k_B * dis;
        dcmplx sum((-1. - I * k_B * dis) * term1 +
                   (-1. + I * k_B * dis) * term2 -
                   dis * 4 * e *
                       exp(-dis * dis * e * e + k_B * k_B / (4 * e * e)) / spi);
        out.second +=
            (sum + 2. * id) * dif / (dis * dis * dis) + k_B * k_B * dif / dis;
      } else {  // Takes limit if exactly zero
        out.first += 2. * I * k_B * (EvErfc(I * k_B / (2 * e)) - id) -
                     4. * e / sqrt(PI) * exp(k_B * k_B / (4. * e * e));
        cvec radial(1., 1., 1.);
        out.second += k_B * k_B * radial / sqrt(3.);
      }
    }
  }
  out.first /= 8 * PI;
  out.second /= -8 * PI;
  return out;
}

// ---------------------------------------------------------------------------
// Ewald split — spectral (reciprocal‑lattice sum)
// ---------------------------------------------------------------------------
// Spectral (reciprocal‑lattice) part of Ewald sum for scalar kernel.
dcmplx GreenFHom3DPer2D::EwaldSpectral(rvec r) {
  dcmplx out(0.);
  rvec rtan(dot(r, a1) / dot(a1, a1) * a1 + dot(r, a2) / dot(a2, a2) * a2);
  rvec rperp(r - rtan);
  double perdif(dot(rperp, perp));
  if (perdif * perdif * e * e < 900.) {
    for (int i(-evaluationTerms); i <= evaluationTerms; ++i) {
      for (int j(-evaluationTerms); j <= evaluationTerms; ++j) {
        rvec kij(k - i * b1 - j * b2);
        dcmplx gamma(csqrt(dot(kij, kij) - k_B * k_B));
        if (imag(gamma) > 0) gamma = -gamma;
        dcmplx gamma2e(gamma / (2 * e));
        dcmplx sum(exp(gamma * perdif) * EvErfc(gamma2e + perdif * e) +
                   exp(-gamma * perdif) * EvErfc(gamma2e - perdif * e));
        out += exp(I * dot(kij, r)) * sum / gamma;
      }
    }
  } else {
    if (perdif < 0) perdif = -perdif;
    for (int i(-evaluationTerms); i <= evaluationTerms; ++i) {
      for (int j(-evaluationTerms); j <= evaluationTerms; ++j) {
        rvec kij(k - i * b1 - j * b2);
        dcmplx gamma(csqrt(dot(kij, kij) - k_B * k_B));
        if (imag(gamma) > 0) gamma = -gamma;
        dcmplx gamma2e(gamma / (2 * e));
        dcmplx sum(exp(-perdif * e * perdif * e) / (sqrt(PI) * perdif * e) +
                   exp(-gamma * perdif) * EvErfc(gamma2e - perdif * e));
        out += exp(I * dot(kij, r)) * sum / gamma;
      }
    }
  }
  out /= 4 * V;
  return out;
}

// Spectral part of Ewald sum for source‑gradient.
cvec GreenFHom3DPer2D::GradientEwaldSpectral(rvec r) {
  cvec out(0.);
  rvec rtan(dot(r, a1) / dot(a1, a1) * a1 + dot(r, a2) / dot(a2, a2) * a2);
  rvec rperp(r - rtan);
  double perdif(dot(rperp, perp));
  if (perdif * perdif * e * e < 900.) {
    for (int i(-evaluationTerms); i <= evaluationTerms; ++i) {
      for (int j(-evaluationTerms); j <= evaluationTerms; ++j) {
        rvec kij(k - i * b1 - j * b2);
        dcmplx gamma(csqrt(dot(kij, kij) - k_B * k_B));
        if (imag(gamma) > 0) gamma = -gamma;
        dcmplx gamma2e(gamma / (2 * e));
        dcmplx sumplus(gamma2e + perdif * e);
        dcmplx summinus(gamma2e - perdif * e);
        dcmplx EvErfcsumplus(EvErfc(sumplus));
        dcmplx EvErfcsumminus(EvErfc(summinus));
        dcmplx expgammaplus(exp(gamma * perdif));
        dcmplx expgammaminus(exp(-gamma * perdif));
        dcmplx phase(exp(I * dot(kij, r)));
        dcmplx sum(expgammaplus * EvErfcsumplus +
                   expgammaminus * EvErfcsumminus);
        out += phase * sum * I * kij / gamma;
        dcmplx secondsum(expgammaplus * EvErfcsumplus -
                         expgammaminus * EvErfcsumminus);
        out += perp * phase * secondsum;
      }
    }
  } else {
    double sign(1.);
    if (perdif < 0) {
      perdif = -perdif;
      sign = -1.;
    }
    for (int i(-evaluationTerms); i <= evaluationTerms; ++i) {
      for (int j(-evaluationTerms); j <= evaluationTerms; ++j) {
        rvec kij(k - i * b1 - j * b2);
        dcmplx gamma(csqrt(dot(kij, kij) - k_B * k_B));
        if (imag(gamma) > 0) gamma = -gamma;
        dcmplx gamma2e(gamma / (2 * e));
        dcmplx summinus(gamma2e - perdif * e);
        dcmplx EvErfcsumminus(EvErfc(summinus));
        dcmplx expgammaminus(exp(-gamma * perdif));
        dcmplx phase(exp(I * dot(kij, r)));
        dcmplx sum(exp(-perdif * e * perdif * e) / (sqrt(PI) * perdif * e) +
                   expgammaminus * EvErfcsumminus);
        out += phase * sum * I * kij / gamma;
        dcmplx secondsum(exp(-perdif * e * perdif * e) /
                             (sqrt(PI) * perdif * e) -
                         expgammaminus * EvErfcsumminus);
        out += sign * perp * phase * secondsum;
      }
    }
  }
  out /= -4 * V;
  return out;
}

// Spectral part returning {G, grad'G} at once.
std::pair<dcmplx, cvec> GreenFHom3DPer2D::EwaldAndGradientSpectral(rvec r) {
  using std::pair;
  pair<dcmplx, cvec> out(0., cvec(0, 0, 0));
  rvec rtan(dot(r, a1) / dot(a1, a1) * a1 + dot(r, a2) / dot(a2, a2) * a2);
  rvec rperp(r - rtan);
  double perdif(dot(rperp, perp));
  if (perdif * perdif * e * e < 900.) {
    for (int i(-evaluationTerms); i <= evaluationTerms; ++i) {
      for (int j(-evaluationTerms); j <= evaluationTerms; ++j) {
        rvec kij(k - i * b1 - j * b2);
        dcmplx gamma(csqrt(dot(kij, kij) - k_B * k_B));
        if (imag(gamma) > 0) gamma = -gamma;
        dcmplx gamma2e(gamma / (2 * e));
        dcmplx sumplus(gamma2e + perdif * e);
        dcmplx summinus(gamma2e - perdif * e);
        dcmplx EvErfcsumplus(EvErfc(sumplus));
        dcmplx EvErfcsumminus(EvErfc(summinus));
        dcmplx expgammaplus(exp(gamma * perdif));
        dcmplx expgammaminus(exp(-gamma * perdif));
        dcmplx phase(exp(I * dot(kij, r)));
        dcmplx sum(expgammaplus * EvErfcsumplus +
                   expgammaminus * EvErfcsumminus);
        out.first += phase * sum / gamma;
        out.second += phase * sum * I * kij / gamma;
        dcmplx secondsum(expgammaplus * EvErfcsumplus -
                         expgammaminus * EvErfcsumminus);
        out.second += perp * phase * secondsum;
      }
    }
  } else {
    double sign(1.);
    if (perdif < 0) {
      perdif = -perdif;
      sign = -1.;
    }
    for (int i(-evaluationTerms); i <= evaluationTerms; ++i) {
      for (int j(-evaluationTerms); j <= evaluationTerms; ++j) {
        rvec kij(k - i * b1 - j * b2);
        dcmplx gamma(csqrt(dot(kij, kij) - k_B * k_B));
        if (imag(gamma) > 0) gamma = -gamma;
        dcmplx gamma2e(gamma / (2 * e));
        dcmplx summinus(gamma2e - perdif * e);
        dcmplx EvErfcsumminus(EvErfc(summinus));
        dcmplx expgammaminus(exp(-gamma * perdif));
        dcmplx phase(exp(I * dot(kij, r)));
        dcmplx sum(exp(-perdif * e * perdif * e) / (sqrt(PI) * perdif * e) +
                   expgammaminus * EvErfcsumminus);
        out.first += phase * sum / gamma;
        out.second += phase * sum * I * kij / gamma;
        dcmplx secondsum(exp(-perdif * e * perdif * e) /
                             (sqrt(PI) * perdif * e) -
                         expgammaminus * EvErfcsumminus);
        out.second += sign * perp * phase * secondsum;
      }
    }
  }
  out.first /= 4 * V;
  out.second /= -4 * V;
  return out;
}

// ---------------------------------------------------------------------------
// Ewald acceleration & pointwise evaluations (Bloch‑periodic)
// ---------------------------------------------------------------------------
/**
 * @details Applies primitive‑cell mapping to the in‑plane displacement and sums
 * spatial & spectral Ewald contributions.
 */
dcmplx GreenFHom3DPer2D::Evaluate(rvec r, rvec rp) { return Ewald(r - rp); }

// Scalar kernel at a fixed lattice translation t.
dcmplx GreenFHom3DPer2D::EvaluateTranslate(rvec r, rvec rp, rvec t) {
  return (EwaldSpatial(r - rp - t) + EwaldSpectral(r - rp - t)) *
         exp(I * dot(k, t));
}

// Ewald sum for the scalar kernel.
dcmplx GreenFHom3DPer2D::Ewald(rvec r) {
  rvec t(PrimitiveCell(r));
  return (EwaldSpatial(r) + EwaldSpectral(r)) * exp(I * dot(k, t));
}

// Smoothed scalar kernel (near singularity) with primitive‑cell mapping.
dcmplx GreenFHom3DPer2D::Smoothed(rvec r, rvec rp) {
  return EwaldSmoothed(r - rp);
}

// Smoothed scalar kernel at fixed translation t.
dcmplx GreenFHom3DPer2D::SmoothedTranslate(rvec r, rvec rp, rvec t) {
  return (EwaldSpatialSmoothed(r - rp - t) + EwaldSpectral(r - rp - t)) *
         exp(I * dot(k, t));
}

// Ewald sum for the smoothed scalar kernel.
dcmplx GreenFHom3DPer2D::EwaldSmoothed(rvec r) {
  rvec t(PrimitiveCell(r));
  return (EwaldSpatialSmoothed(r) + EwaldSpectral(r)) * exp(I * dot(k, t));
}

// Source‑gradient grad'G for periodic kernel with primitive‑cell mapping.
cvec GreenFHom3DPer2D::Gradient(rvec r, rvec rp) {
  return GradientEwald(r - rp);
}

// Source‑gradient grad'G at fixed translation t.
cvec GreenFHom3DPer2D::GradientTranslate(rvec r, rvec rp, rvec t) {
  return (GradientEwaldSpatial(r - rp - t) +
          GradientEwaldSpectral(r - rp - t)) *
         exp(I * dot(k, t));
}

// Ewald sum for the gradient w.r.t. source (grad'G).
cvec GreenFHom3DPer2D::GradientEwald(rvec r) {
  rvec t(PrimitiveCell(r));
  return (GradientEwaldSpatial(r) + GradientEwaldSpectral(r)) *
         exp(I * dot(k, t));
}

// Smoothed source‑gradient with primitive‑cell mapping.
cvec GreenFHom3DPer2D::GradientSmoothed(rvec r, rvec rp) {
  return GradientEwaldSmoothed(r - rp);
}

// Smoothed source‑gradient at fixed translation t.
cvec GreenFHom3DPer2D::GradientSmoothedTranslate(rvec r, rvec rp, rvec t) {
  return (GradientEwaldSpatialSmoothed(r - rp - t) +
          GradientEwaldSpectral(r - rp - t)) *
         exp(I * dot(k, t));
}

// Ewald sum for the smoothed gradient w.r.t. source (grad'G). 
cvec GreenFHom3DPer2D::GradientEwaldSmoothed(rvec r) {
  rvec t(PrimitiveCell(r));
  return (GradientEwaldSpatialSmoothed(r) + GradientEwaldSpectral(r)) *
         exp(I * dot(k, t));
}

// Evaluate {G, grad'G} together at a fixed lattice translation t.
std::pair<dcmplx, cvec> GreenFHom3DPer2D::EvaluateandGradientTranslate(rvec r,
                                                                       rvec rp,
                                                                       rvec t) {
  std::pair<dcmplx, cvec> out1, out2, out;
  out1 = EwaldAndGradientSpatial(r - rp - t);
  out2 = EwaldAndGradientSpectral(r - rp - t);
  out.first = (out1.first + out2.first) * exp(I * dot(k, t));
  out.second = (out1.second + out2.second) * exp(I * dot(k, t));
  return out;
}

// Evaluate smoothed {G, grad'G} together at fixed translation t.
std::pair<dcmplx, cvec> GreenFHom3DPer2D::EvaluateandGradientSmoothedTranslate(
    rvec r, rvec rp, rvec t) {
  std::pair<dcmplx, cvec> out1, out2, out;
  out1 = EwaldAndGradientSpatialSmoothed(r - rp - t);
  out2 = EwaldAndGradientSpectral(r - rp - t);
  out.first = (out1.first + out2.first) * exp(I * dot(k, t));
  out.second = (out1.second + out2.second) * exp(I * dot(k, t));
  return out;
}

// ---------------------------------------------------------------------------
// RWG‑weighted integrals
// ---------------------------------------------------------------------------
// Triangle–triangle D entry (disabled for periodic homogeneous case).
dcmplx GreenFHom3DPer2D::IntegrateD(RWGFun* /*f*/, RWGFun* /*fp*/) { exit(1); }

// Triangle–triangle K entry (disabled for periodic homogeneous case).
dcmplx GreenFHom3DPer2D::IntegrateK(RWGFun* /*f*/, RWGFun* /*fp*/) { exit(1); }

/**
 * @details Uses Evaluate/Smoothed and Gradient/GradientSmoothed according to the
 * near/far heuristic, and applies analytic singularity subtraction.
 */
cvec GreenFHom3DPer2D::IntegrateDelta(rvec r, RWGFun* fp) {
  cvec out(0., 0., 0.);
  Triangle* Tp = fp->TrianglePtr();
  rvec rp[3];
  for (int i = 0; i < 3; i++) {
    rp[i] = Tp->Node(i);
  }
  double Ap = Tp->Area(), perimp = Tp->Perimeter();
  rvec Tc(Tp->Center());
  rvec junk(r - Tc);
  rvec t(PrimitiveCell(junk));
  double dist = sqrt(dot(r - t - Tc, r - t - Tc));
  double kr = sqrt(norm(k_B)) * dist;

  int Ndunavant;
  double(*xdunavant)[3], *wdunavant;
  if (NeedsAccurate) {
    Ndunavant = (dist < perimp) ? dunavant::N17 : dunavant::N5;
    xdunavant = (dist < perimp) ? dunavant::x17 : dunavant::x5;
    wdunavant = (dist < perimp) ? dunavant::w17 : dunavant::w5;
  } else {
    Ndunavant = dunavant::N1;
    xdunavant = dunavant::x1;
    wdunavant = dunavant::w1;
  }

  bool NeedsSingSub = (dist < perimp) || (kr < 0.2 * PI);
  std::pair<dcmplx, cvec> (GreenFHom3DPer2D::*DKfunction)(rvec r, rvec rp,
                                                          rvec t) =
      NeedsSingSub ? &GreenFHom3DPer2D::EvaluateandGradientSmoothedTranslate
                   : &GreenFHom3DPer2D::EvaluateandGradientTranslate;

  std::pair<dcmplx, cvec> DKtemp;
  dcmplx Dtemp;
  cvec Ktemp;
  for (int i = 0; i < Ndunavant; ++i) {
    rvec sp = rp[0] * xdunavant[i][0] + rp[1] * xdunavant[i][1] +
              rp[2] * xdunavant[i][2];
    double dAp = Ap * wdunavant[i];
    DKtemp = (this->*DKfunction)(r, sp, t);
    Dtemp = dAp * DKtemp.first;
    Ktemp = dAp * DKtemp.second;
    rvec fpsp(fp->Evaluate(sp));
    out += (fpsp * Dtemp - 1. / (k_B * k_B) * fp->RWGPreFactor() * Ktemp);
  }

  if (NeedsSingSub) {
    SingSub ss;
    ss.SetTriangleRWG(fp, r - t);
    out += exp(I * dot(this->k, t)) / (4. * PI) *
           (-1. / (k_B * k_B) * fp->RWGPreFactor() *
                (ss.K3(-1) - k_B * k_B / 2. * ss.K3(1)) +
            (ss.K2(-1) - k_B * k_B / 2. * ss.K2(1)));
  }

  return out;
}

// Point–triangle Κ vector integral at observation r (2D‑periodic kernel).
cvec GreenFHom3DPer2D::IntegrateKappa(rvec r, RWGFun* fp) {
  cvec out(0., 0., 0.);
  Triangle* Tp = fp->TrianglePtr();
  rvec rp[3];
  for (int i = 0; i < 3; i++) {
    rp[i] = Tp->Node(i);
  }
  double Ap = Tp->Area(), perimp = Tp->Perimeter();
  rvec Tc(Tp->Center());
  rvec junk(r - Tc);
  rvec t(PrimitiveCell(junk));
  double dist = sqrt(dot(r - t - Tc, r - t - Tc));
  double kr = sqrt(norm(k_B)) * dist;

  int Ndunavant;
  double(*xdunavant)[3], *wdunavant;
  if (NeedsAccurate) {
    Ndunavant = (dist < perimp) ? dunavant::N17 : dunavant::N5;
    xdunavant = (dist < perimp) ? dunavant::x17 : dunavant::x5;
    wdunavant = (dist < perimp) ? dunavant::w17 : dunavant::w5;
  } else {
    Ndunavant = dunavant::N1;
    xdunavant = dunavant::x1;
    wdunavant = dunavant::w1;
  }

  bool NeedsSingSub = (dist < perimp) || (kr < 0.2 * PI);
  cvec (GreenFHom3DPer2D::*Kfunction)(rvec r, rvec rp, rvec t) =
      NeedsSingSub ? &GreenFHom3DPer2D::GradientSmoothedTranslate
                   : &GreenFHom3DPer2D::GradientTranslate;

  cvec Ktemp;
  for (int i = 0; i < Ndunavant; ++i) {
    rvec sp = rp[0] * xdunavant[i][0] + rp[1] * xdunavant[i][1] +
              rp[2] * xdunavant[i][2];
    double dAp = Ap * wdunavant[i];
    Ktemp = dAp * (this->*Kfunction)(r, sp, t);
    rvec fpsp(fp->Evaluate(sp));
    out += cross(Ktemp, (cvec)fpsp);
  }

  if (NeedsSingSub) {
    SingSub ss;
    ss.SetTriangleRWG(fp, r - t);
    out += exp(I * dot(this->k, t)) / (4. * PI) *
           (ss.K4(-1) - k_B * k_B / 2. * ss.K4(1));
  }

  return out;
}

std::pair<dcmplx, dcmplx> GreenFHom3DPer2D::PairD(RWGFun* f, RWGFun* fp) {
  // Don't care about symmetry, this is more exact!
  dcmplx elem(IntegrateD(f, fp));
  if (f == fp) return std::pair<dcmplx, dcmplx>(elem, elem);
  dcmplx elem2(IntegrateD(fp, f));
  return std::pair<dcmplx, dcmplx>(elem, elem2);
}

std::pair<dcmplx, dcmplx> GreenFHom3DPer2D::PairK(RWGFun* f, RWGFun* fp) {
  // Unfortunately, IntegrateK is not symmetric.
  dcmplx elem(IntegrateK(f, fp));
  if (f == fp) return std::pair<dcmplx, dcmplx>(elem, elem);
  dcmplx elem2(IntegrateK(fp, f));
  return std::pair<dcmplx, dcmplx>(elem, elem2);
}

// Batched D/K block assembly for same‑triangle RWG sets (periodic case).
void GreenFHom3DPer2D::SameTriDK(
    const std::vector<RWGFun*>& fvec, const std::vector<RWGFun*>& fpvec,
    std::vector<std::vector<std::pair<dcmplx, dcmplx> > >& DK) {
  Triangle *tPtr = fvec[0]->TrianglePtr(), *tpPtr = fpvec[0]->TrianglePtr();
  rvec r[3], rp[3];
  for (int i = 0; i < 3; i++) {
    r[i] = tPtr->Node(i);
    rp[i] = tpPtr->Node(i);
  }
  double A = tPtr->Area(), Ap = tpPtr->Area();
  SingSub ss;

  // Primitive cell is being found only once based on the triangle centres
  // All other functions are modified to pass this directly
  rvec junk(tPtr->Center() - tpPtr->Center());
  rvec t(PrimitiveCell(junk));
  int tristatus = tPtr->FindAdjacency(tpPtr, t);

  int Ndunavant;
  double(*xdunavant)[3];
  double* wdunavant;
  bool NeedsDanalytic(false);
  bool NeedsDSingSub(false);
  bool NeedsKSingSub(false);
  bool NeedsKLineIntegral(false);
  bool DistantAboveThreshold;
  std::pair<dcmplx, cvec> (GreenFHom3DPer2D::*DKfunction)(rvec r, rvec rp,
                                                          rvec t);
  if (NeedsAccurate) {
    Ndunavant = tristatus > 0 ? dunavant::N17 : dunavant::N5;
    xdunavant = tristatus > 0 ? dunavant::x17 : dunavant::x5;
    wdunavant = tristatus > 0 ? dunavant::w17 : dunavant::w5;
    NeedsDanalytic = tristatus == 3;
    NeedsKLineIntegral = tristatus == 2;
    DistantAboveThreshold =
        tristatus == 0 and AboveThresholdTranslate(tPtr, tpPtr, t);
  } else {
    Ndunavant = dunavant::N1;
    xdunavant = dunavant::x1;
    wdunavant = dunavant::w1;
    DistantAboveThreshold = AboveThreshold(tPtr->Center() - t, tpPtr->Center());
  }
  NeedsDSingSub = !DistantAboveThreshold;
  NeedsKSingSub = tristatus != 3 and !DistantAboveThreshold;
  DKfunction = DistantAboveThreshold
                   ? &GreenFHom3DPer2D::EvaluateandGradientTranslate
                   : &GreenFHom3DPer2D::EvaluateandGradientSmoothedTranslate;

  std::vector<rvec> svec(Ndunavant), spvec(Ndunavant);
  for (int i = 0; i < Ndunavant; ++i) {
    svec[i] = r[0] * xdunavant[i][0] + r[1] * xdunavant[i][1] +
              r[2] * xdunavant[i][2];
    spvec[i] = rp[0] * xdunavant[i][0] + rp[1] * xdunavant[i][1] +
               rp[2] * xdunavant[i][2];
  }

  std::pair<dcmplx, cvec> DKtemp;
  dcmplx Dtemp;
  cvec Ktemp;
  for (int i = 0; i < Ndunavant; ++i) {
    double dA = A * wdunavant[i];
    for (int j = 0; j < Ndunavant; ++j) {
      double dAp = Ap * wdunavant[j];
      DKtemp = (this->*DKfunction)(svec[i], spvec[j], t);
      Dtemp = dA * dAp * DKtemp.first;
      Ktemp = dA * dAp * DKtemp.second;
      for (int k = 0; k < (int)fvec.size(); ++k) {
        RWGFun* f = fvec[k];
        rvec fs(f->Evaluate(svec[i]));
        double divf(f->RWGPreFactor());
        for (int l = 0; l < (int)fpvec.size(); ++l) {
          RWGFun* fp = fpvec[l];
          rvec fpsp(fp->Evaluate(spvec[j]));
          double divfp(fp->RWGPreFactor());

          DK[k][l].first +=
              Dtemp * (dot(fs, fpsp) - 1. / (k_B * k_B) * divf * divfp);
          DK[k][l].second += dot(fs, cross(Ktemp, (cvec)fpsp));
        }
      }
    }

    if (NeedsDSingSub or NeedsKSingSub) {
      ss.SetTriangleRWG(fpvec[0], svec[i] - t);
      for (int l = 0; l < (int)fpvec.size(); ++l) {
        RWGFun* fp = fpvec[l];
        double divfp(fp->RWGPreFactor());
        ss.ChangeTriangleRWG(fp);
        for (int k = 0; k < (int)fvec.size(); ++k) {
          RWGFun* f = fvec[k];
          double divf(f->RWGPreFactor());
          rvec fs(f->Evaluate(svec[i]));

          if (NeedsDanalytic) {
            DK[k][l].first +=
                exp(I * dot(this->k, t)) * dA / (4. * PI) *
                (-1. / (k_B * k_B) * divf * (-k_B * k_B / 2. * ss.K1(1)) +
                 dot(fs, -k_B * k_B / 2. * ss.K2(1)));
          } else if (NeedsDSingSub) {
            DK[k][l].first += exp(I * dot(this->k, t)) * dA / (4. * PI) *
                              (-1. / (k_B * k_B) * divf *
                                   (ss.K1(-1) - k_B * k_B / 2. * ss.K1(1)) +
                               dot(fs, ss.K2(-1) - k_B * k_B / 2. * ss.K2(1)));
          }
          if (NeedsKLineIntegral) {
            DK[k][l].second +=
                exp(I * dot(this->k, t)) * dA / (4. * PI) *
                dot(fs, cross((rvec)(f->FreeVertex() - fp->FreeVertex() - t),
                              ss.K3n(-1)) *
                                divfp / (-2.) -
                            k_B * k_B / 2. * ss.K4(1));
          } else if (NeedsKSingSub) {
            DK[k][l].second += exp(I * dot(this->k, t)) * dA / (4. * PI) *
                               dot(fs, ss.K4(-1) - k_B * k_B / 2. * ss.K4(1));
          }
        }
      }
    }
  }

  if (NeedsKLineIntegral) {
    // Line integral for surface gradient of 1/R
    // Note that triangle order is reversed
    rvec n = cross((rvec)(rp[1] - rp[0]), (rvec)(rp[2] - rp[0]));
    for (int i = 0; i < 3; i++) {
      rvec &p1 = rp[i], &p2 = rp[(i + 1) % 3];
      double L = sqrt(dot(p2 - p1, p2 - p1));

      // Outward normal on edge;
      rvec m = cross((rvec)(p2 - p1), n);
      m = m / sqrt(dot(m, m));

      // Line integral
      for (int j = 0; j < gausslegendre::N17; j++) {
        rvec s = p1 * gausslegendre::x17[j][0] + p2 * gausslegendre::x17[j][1];
        double dL = L * gausslegendre::w17[j];

        ss.SetTriangleRWG(fvec[0], s + t);
        for (int k = 0; k < (int)fvec.size(); ++k) {
          RWGFun* f = fvec[k];
          ss.ChangeTriangleRWG(f);
          for (int l = 0; l < (int)fpvec.size(); ++l) {
            RWGFun* fp = fpvec[l];
            rvec pqm = cross((rvec)(f->FreeVertex() - fp->FreeVertex() - t), m);

            DK[k][l].second += -exp(I * dot(this->k, t)) * dL / (4. * PI) *
                               fp->RWGPreFactor() / 2. * dot(pqm, ss.K2(-1));
          }
        }
      }
    }
  }

  if (NeedsDanalytic) {
    for (int k = 0; k < (int)fvec.size(); ++k) {
      RWGFun* f = fvec[k];
      double divf(f->RWGPreFactor());
      for (int l = 0; l < (int)fpvec.size(); ++l) {
        RWGFun* fp = fpvec[l];
        double divfp(fp->RWGPreFactor());
        // Analytical 4D integral for (div f)(div f')/(4 pi r)
        DK[k][l].first += -exp(I * dot(this->k, t)) * 1. / (k_B * k_B) * divf *
                          divfp / (4 * PI) * tPtr->getI1();
        // Analytical 4D integral for f dot f'/(4 pi r)
        // Note that f=fp and f!=fp need to be treated differently
        DK[k][l].first +=
            exp(I * dot(this->k, t)) * divf * divfp / (16. * PI) *
            tpPtr->getI2Translated(f->FreeVertex(), fp->FreeVertex(), t);
      }
    }
  }
}

// Batched delta/kappa evaluation for RWGs on one triangle (periodic case).
void GreenFHom3DPer2D::SameTriDeltaKappa(
    rvec r, const std::vector<RWGFun*>& fpvec,
    std::vector<std::pair<cvec, cvec> >& DeltaKappa) {
  Triangle* Tp = fpvec[0]->TrianglePtr();
  rvec rp[3];
  for (int i = 0; i < 3; i++) {
    rp[i] = Tp->Node(i);
  }
  double Ap = Tp->Area(), perimp = Tp->Perimeter();
  rvec Tc(Tp->Center());
  rvec junk(r - Tc);
  rvec t(PrimitiveCell(junk));
  double dist = sqrt(dot(r - t - Tc, r - t - Tc));
  double kr = sqrt(norm(k_B)) * dist;
  SingSub ss;

  int Ndunavant;
  double(*xdunavant)[3], *wdunavant;
  if (NeedsAccurate) {
    Ndunavant = (dist < perimp) ? dunavant::N17 : dunavant::N5;
    xdunavant = (dist < perimp) ? dunavant::x17 : dunavant::x5;
    wdunavant = (dist < perimp) ? dunavant::w17 : dunavant::w5;
  } else {
    Ndunavant = dunavant::N1;
    xdunavant = dunavant::x1;
    wdunavant = dunavant::w1;
  }

  bool NeedsSingSub = (dist < perimp) || (kr < 0.2 * PI);
  std::pair<dcmplx, cvec> (GreenFHom3DPer2D::*DKfunction)(rvec r, rvec rp,
                                                          rvec t) =
      NeedsSingSub ? &GreenFHom3DPer2D::EvaluateandGradientSmoothedTranslate
                   : &GreenFHom3DPer2D::EvaluateandGradientTranslate;

  std::pair<dcmplx, cvec> DKtemp;
  dcmplx Dtemp;
  cvec Ktemp;
  for (int i = 0; i < Ndunavant; ++i) {
    rvec sp = rp[0] * xdunavant[i][0] + rp[1] * xdunavant[i][1] +
              rp[2] * xdunavant[i][2];
    double dAp = Ap * wdunavant[i];
    DKtemp = (this->*DKfunction)(r, sp, t);
    Dtemp = dAp * DKtemp.first;
    Ktemp = dAp * DKtemp.second;
    for (int j = 0; j < (int)fpvec.size(); ++j) {
      RWGFun* fp = fpvec[j];
      rvec fpsp(fp->Evaluate(sp));
      DeltaKappa[j].first +=
          (fpsp * Dtemp - 1. / (k_B * k_B) * fp->RWGPreFactor() * Ktemp);
      DeltaKappa[j].second += cross(Ktemp, (cvec)fpsp);
    }
  }

  if (NeedsSingSub) {
    ss.SetTriangleRWG(fpvec[0], r - t);
    for (int j = 0; j < (int)fpvec.size(); ++j) {
      RWGFun* fp = fpvec[j];
      ss.ChangeTriangleRWG(fp);
      DeltaKappa[j].first += exp(I * dot(this->k, t)) / (4. * PI) *
                             (-1. / (k_B * k_B) * fp->RWGPreFactor() *
                                  (ss.K3(-1) - k_B * k_B / 2. * ss.K3(1)) +
                              (ss.K2(-1) - k_B * k_B / 2. * ss.K2(1)));
      DeltaKappa[j].second += exp(I * dot(this->k, t)) / (4. * PI) *
                              (ss.K4(-1) - k_B * k_B / 2. * ss.K4(1));
    }
  }
}
