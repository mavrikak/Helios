/**
 * @file Dipole.cpp
 * @brief Implementation of the #Dipole incident source.
 */

#include "Dipole.h"
#include <iostream>

// ---------------------------
// Constructors
// ---------------------------
Dipole::Dipole(rvec location, cvec polarization, dcmplx vacuumWavelength,
               dcmplx mu, GreenF* fHom)
    : p(polarization),
      r0(location),
      omega(2 * PI * CVAC / vacuumWavelength),
      grnFun(fHom),
      muMedium(mu) { grnFunLayered = nullptr; grnFunLayeredFields = nullptr;}

Dipole::Dipole(rvec location, cvec polarization, dcmplx vacuumWavelength, 
               GreenF* fHom, dcmplx mu, std::vector<Grid> tabGrids, 
               GreenF* fLayered, LayeredMediaUtils* utils)
    : p(polarization),
      r0(location),
      omega(2 * PI * CVAC / vacuumWavelength),
      grnFun(fHom), 
      muMedium(mu),
      tabulationGrids(tabGrids),
      grnFunLayered(fLayered),
      layeredUtils(utils) { grnFunLayeredFields = nullptr; }

// ---------------------------
// Pointwise fields (homogeneous)
// ---------------------------
cvec Dipole::EvaluateE(rvec r) {
  return omega * omega * MU0 * muMedium * dot(grnFun->EvaluateDyadic(r, r0), p);
}

cvec Dipole::EvaluateH(rvec r) {
  // Minus because Gradient() is with respect to second parameter
  return I * omega * cross(grnFun->Gradient(r, r0), p);
}

// ---------------------------
// Pointwise fields (layered)
// ---------------------------
std::vector<cvec> Dipole::EvaluateLayeredEandH(rvec r)
{
  // Start with zero fields
  cvec Einc(dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0));
  cvec Hinc(dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0));

  // Homogeneous contribution added only if r and r0 reside in the same layer
  int obsLayerIndex = this->layeredUtils->GetLayerIndex(r);
  int incLayerIndex = this->layeredUtils->GetLayerIndex(r0);
  if (obsLayerIndex == incLayerIndex) {
    Einc += this->EvaluateE(r);
    Hinc += this->EvaluateH(r);
  }
  
  // Layered Green's function contributions
  if (!grnFunLayeredFields) {
    throw std::runtime_error("Green function not set on this Dipole.");
  }
  // Downcast to GreenFLayered3D to access specific methods
  auto* grnPtr = dynamic_cast<GreenFLayered3D*>(grnFunLayeredFields);
  if (!grnPtr) {
    throw std::runtime_error("Green function pointer is not a GreenFLayered3D.");
  }

  // Set dyadic type to "normal" for E field contribution
  Einc += omega * omega * MU0 * muMedium * dot(grnPtr->EvaluateDyadic(r, r0, "normal"), p);

  // Set dyadic type to "curl" for H field contribution
  Hinc += I * omega * dot(grnPtr->EvaluateDyadic(r, r0, "curl"), p);
  
  return {Einc, Hinc};
}

// ---------------------------
// RWG integrations (homogeneous and layered)
// ---------------------------
dcmplx Dipole::IntegrateE(RWGFun* f) {
  return omega * omega * MU0 * muMedium * dot(p, grnFun->IntegrateDelta(r0, f));
}

dcmplx Dipole::IntegrateH(RWGFun* f) {
  return I * omega * dot(p, grnFun->IntegrateKappa(r0, f));
}

// Combine direct (homogeneous) and reflected/transmitted (layered) parts.
std::vector<dcmplx> Dipole::IntegrateLayeredEandH(RWGFun* f) {
  dcmplx Einc = omega * omega * MU0 * muMedium *
                ( dot( p, grnFun->IntegrateDelta( r0, f ) ) + 
                  dot( p, grnFunLayered->IntegrateDelta( r0, f ) ) );
  dcmplx Hinc = I * omega * 
                ( dot( p, grnFun->IntegrateKappa( r0, f ) ) + 
                  dot( p, grnFunLayered->IntegrateKappa( r0, f ) ) );
  return {Einc, Hinc};
}

// ---------------------------
// Accessors & configuration
// ---------------------------
rvec Dipole::Location() { return r0; }

cvec Dipole::Polarization() { return p; }

// Set Green's function for layered media
void Dipole::setGrnFunLayered(GreenF* fGrnFields, std::vector<Grid> tabGrids) {
  grnFunLayeredFields = fGrnFields;
  tabulationGridsFields = tabGrids;
}