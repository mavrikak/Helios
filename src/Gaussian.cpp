/**
 * @file Gaussian.cpp
 * @brief Implementation of the paraxial Gaussian beam incident field.
 */

#include "Gaussian.h"
#include <cmath>
#include <iostream>

using namespace std;

// -----------------------------
// Construction
// -----------------------------
Gaussian::Gaussian(double wavelength, double waist, rvec focal_point,
                   double eps, rvec propdir, cvec polarisation) {
  prop = propdir / sqrt(dot(propdir, propdir));
  k = (2 * PI * sqrt(eps) / wavelength);
  w0 = waist;
  focus = focal_point;
  z = Z0 * (sqrt(eps));
  p = polarisation;
  depth = (PI * w0 * w0 * sqrt(eps)) / wavelength;
}

// -----------------------------
// Accessors
// -----------------------------
double Gaussian::Waist() const { return w0; }

cvec Gaussian::Polarisation() const { return p; }

rvec Gaussian::Propagation() const { return k * prop; }

// -----------------------------
// Pointwise field evaluation
// -----------------------------
cvec Gaussian::EvaluateE(rvec r) {
  // Centering the coordinate system in the focus point
  r = r - focus;

  // distance from focus along the propagation axis
  double Z = dot(r, prop);

  // distance from z axis
  double rho = sqrt(dot(r - Z * prop, r - Z * prop));

  // beam size at position r
  double w = w0 * sqrt(1 + (Z / depth) * (Z / depth));

  // curvature of wavefront at position r
  double C = Z / (Z * Z + depth * depth);

  // Gouy phase
  double zeta = atan(Z / depth);

  // Electric field of a paraxial gaussian beam
  return p * (w0 / w) * exp(-rho * rho / (w * w) +
                            I * k * (Z + C * rho * rho * 0.5) - I * zeta);
}

cvec Gaussian::EvaluateH(rvec r) {
  cvec propagation(prop);
  return (1 / z) * cross(propagation, EvaluateE(r));
}

std::vector<cvec> Gaussian::EvaluateLayeredEandH(rvec /*r*/) {
  throw std::runtime_error("Gaussian::EvaluateLayeredEandH not implemented.");
}

// -----------------------------
// RWG‑weighted integrals
// -----------------------------
dcmplx Gaussian::IntegrateE(RWGFun* f) {
  dcmplx out;
  cvec (Gaussian::*Evaluationfct)(rvec r);
  Evaluationfct = &Gaussian::EvaluateE;
  out = gauss.Integrate(this, Evaluationfct, f);
  return out;
}

dcmplx Gaussian::IntegrateH(RWGFun* f) {
  dcmplx out;
  cvec (Gaussian::*Evaluationfct)(rvec r);
  Evaluationfct = &Gaussian::EvaluateH;
  out = gauss.Integrate(this, Evaluationfct, f);
  return out;
}

std::vector<dcmplx> Gaussian::IntegrateLayeredEandH(RWGFun* /*f*/) {
  throw std::runtime_error("Gaussian::IntegrateLayeredEandH not implemented.");
}

// -----------------------------
// Layered‑media wiring (unsupported)
// -----------------------------
void Gaussian::setGrnFunLayered(GreenF* /*fGrnFields*/, std::vector<Grid> /*tabGrids*/) {
  throw std::runtime_error("Gaussian does not support layered Green's functions.");
}