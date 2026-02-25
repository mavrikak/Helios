/**
 * @file DomainHom3D.cpp
 * @brief Implementation of the homogeneous 3D domain.
 */

#include "DomainHom3D.h"
#include <iostream>
#include "Dipole.h"
#include "IncidentField.h"
#include "PlaneWave.h"
#include "SIEFormPMCHW.h"
#include "iofunctions.h"
#include "Gaussian.h"

// -----------------------------
// Construction & lifetime
// -----------------------------
DomainHom3D::DomainHom3D(SurfaceMesh *inMesh, int inDIndex, dcmplx inEpsilon,
                         dcmplx inMu, dcmplx inVacuumWavelength)
    : Domain(inMesh, inDIndex, inVacuumWavelength),
      epsilon(inEpsilon),
      mu(inMu) {
  InitGrnFun();
}

DomainHom3D::DomainHom3D(SurfaceMesh *inMesh, int inDIndex, dcmplx inEpsilon,
                         dcmplx inMu, dcmplx inVacuumWavelength, bool parent)
    : Domain(inMesh, inDIndex, inVacuumWavelength),
      epsilon(inEpsilon),
      mu(inMu) {
  if (!parent) InitGrnFun();
}

DomainHom3D::~DomainHom3D() { delete grnFun; }

// -----------------------------
// Green's function & formulation
// -----------------------------
int DomainHom3D::InitGrnFun() {
  grnFun = new GreenFHom3D(vacuumWavelength, epsilon, mu);
  return 0;
}

int DomainHom3D::InitFormulation(dcmplx inVacWavelength) {
  formulation = new SIEFormPMCHW(&rwgFuns, grnFun, &incidentFields,
                                 inVacWavelength, epsilon, mu);
  return 0;
}

// -----------------------------
// Cross‑sections
// -----------------------------
/**
 * @brief Compute the element-wise complex conjugate of a 3-vector.
 * @param in Input complex vector.
 * @return Complex vector containing conjugated components of @p in.
 */
cvec conjVecHom(cvec in) { return cvec(conj(in(0)), conj(in(1)), conj(in(2))); }

/**
 * @brief Extract the real part of each component of a complex 3-vector.
 * @param in Input complex vector.
 * @return Real vector with the real parts of @p in.
 */
rvec realVecHom(cvec in) { return rvec(real(in(0)), real(in(1)), real(in(2))); }

/**
 * @details Builds per‑triangle field samples (E,H) from RWG coefficients,
 * integrates the outward Poynting component of the scattered field, and
 * the near-field intensity.
 */
std::vector<double> DomainHom3D::scatteringCS(blitz::Array<dcmplx,1>* solVector)
{
  // Gather unique triangles & normals
  std::vector<Triangle*> vTri;
  std::vector<rvec> vNormals;
  vTri.reserve(rwgFuns.size());
  vNormals.reserve(rwgFuns.size());

  for (RWGFun& rwg : rwgFuns) {
    Triangle* T = rwg.TrianglePtr();
    auto it = std::find(vTri.begin(), vTri.end(), T);
    if (it == vTri.end()) {
      vTri.push_back(T);
      vNormals.push_back(T->FBNormalVector(rwg.Side()));
    }
  }

  size_t Ntri = vTri.size();  // number of unique triangles

  // Quadrature rule on reference triangle
  int nQuad = dunavant::N5;    // number of quadrature points
  auto& xQuad = dunavant::x5;  // barycentrics [nQuad][3]
  auto& wQuad = dunavant::w5;  // weights[nQuad]

  // Allocate per-triangle, per-quad-point accumulators
  std::vector<std::vector<cvec>> Etot(Ntri, std::vector<cvec>(nQuad, cvec(0, 0, 0))),
                                 Htot(Ntri, std::vector<cvec>(nQuad, cvec(0, 0, 0)));

  int offset = solVector->shape()(0) / 2; // offset for β coefficients

  // Fill E, H normal and tangential components at each sample point
  for (RWGFun& rwg : rwgFuns) {
    // find which triangle index
    Triangle* T = rwg.TrianglePtr();
    size_t triIndex = std::distance(vTri.begin(),
                      std::find(vTri.begin(), vTri.end(), T));
    rvec normal = vNormals[triIndex];

    // RWG coefficients
    dcmplx alpha = (*solVector)(rwg.EdgeIndex());
    dcmplx beta  = (*solVector)(rwg.EdgeIndex() + offset);
    double divf  = rwg.RWGPreFactor();

    // triangle nodes
    rvec r[3] = { T->Node(0), T->Node(1), T->Node(2) };

    // loop quad points
    for (int q = 0; q < nQuad; ++q) {
      // physical sample point
      rvec s = xQuad[q][0] * r[0] + xQuad[q][1] * r[1] + xQuad[q][2] * r[2];

      // hRWG at s
      rvec fval = rwg.Evaluate(s);

      // normal‐component (constant over tri)
      cvec En_contrib = -I * alpha * divf * normal / (omega * epsilon * EPS0);
      cvec Hn_contrib = -I * beta * divf * normal / (omega * mu * MU0);

      // tangential‐component
      cvec Et_contrib =  cross( (cvec)normal, (cvec)(beta * fval));
      cvec Ht_contrib = -cross( (cvec)normal, (cvec)(alpha * fval));

      // J,M currents have contributions from each RWG with opposite sign, thus -=
      Etot[triIndex][q] -= En_contrib + Et_contrib;
      Htot[triIndex][q] -= Hn_contrib + Ht_contrib;
    }
  }

  // Quadrature - find the scattering cross section
  double sigma = 0.0;
  double near = 0.0;
  for (size_t t = 0; t < Ntri; ++t) {
    rvec normal = vNormals[t];
    double area = vTri[t]->Area();
    
    // re-compute nodes for this triangle
    rvec r[3] = { vTri[t]->Node(0), vTri[t]->Node(1), vTri[t]->Node(2) };

    for (int q = 0; q < nQuad; ++q) {
      // map barycentrics -> point
      rvec s = xQuad[q][0] * r[0] + xQuad[q][1] * r[1] + xQuad[q][2] * r[2];

      cvec Einc = IncFieldE(s);
      cvec Hinc = IncFieldH(s);
      rvec Sinc = realVecHom( cross(Einc, conjVecHom(Hinc)) ) / 2.0;

      cvec Escat = Etot[t][q] - Einc;
      cvec Hscat = Htot[t][q] - Hinc;
      rvec Sscat = realVecHom( cross(Escat, conjVecHom(Hscat)) ) / 2.0;

      sigma += wQuad[q] * area * dot(normal, Sscat) / sqrt( dot( Sinc, Sinc ) );
      near += real( wQuad[q] * area * dot( Escat, conjVecHom(Escat) ) );
    }
  }
  return {sigma, near};
}

/**
 * @details Integrates the inward Poynting component of total fields.
 */
double DomainHom3D::absorptionCS(blitz::Array<dcmplx,1>* solVector)
{
  // Gather unique triangles & normals
  std::vector<Triangle*> vTri;
  std::vector<rvec> vNormals;
  vTri.reserve(rwgFuns.size());
  vNormals.reserve(rwgFuns.size());

  for (RWGFun& rwg : rwgFuns) {
    Triangle* T = rwg.TrianglePtr();
    auto it = std::find(vTri.begin(), vTri.end(), T);
    if (it == vTri.end()) {
      vTri.push_back(T);
      vNormals.push_back(T->FBNormalVector(rwg.Side()));
    }
  }

  size_t Ntri = vTri.size();  // number of unique triangles

  // Quadrature rule on reference triangle
  int nQuad = dunavant::N5;    // number of quadrature points
  auto& xQuad = dunavant::x5;  // barycentrics [nQuad][3]
  auto& wQuad = dunavant::w5;  // weights[nQuad]

  // Allocate per-triangle, per-quad-point accumulators
  std::vector<std::vector<cvec>> Etot(Ntri, std::vector<cvec>(nQuad, cvec(0, 0, 0))),
                                 Htot(Ntri, std::vector<cvec>(nQuad, cvec(0, 0, 0)));

  int offset = solVector->shape()(0) / 2; // offset for β coefficients

  // Fill E, H normal and tangential components at each sample point
  for (RWGFun& rwg : rwgFuns) {
    // find which triangle index
    Triangle* T = rwg.TrianglePtr();
    size_t triIndex = std::distance(vTri.begin(),
                      std::find(vTri.begin(), vTri.end(), T));
    rvec normal = vNormals[triIndex];

    // RWG coefficients
    dcmplx alpha = (*solVector)(rwg.EdgeIndex());
    dcmplx beta  = (*solVector)(rwg.EdgeIndex() + offset);
    double divf  = rwg.RWGPreFactor();

    // triangle nodes
    rvec r[3] = { T->Node(0), T->Node(1), T->Node(2) };

    // loop quad points
    for (int q = 0; q < nQuad; ++q) {
      // physical sample point
      rvec s = xQuad[q][0] * r[0] + xQuad[q][1] * r[1] + xQuad[q][2] * r[2];

      // hRWG at s
      rvec fval = rwg.Evaluate(s);

      // normal‐component (constant over tri)
      cvec En_contrib = -I * alpha * divf * normal / (omega * epsilon * EPS0);
      cvec Hn_contrib = -I * beta * divf * normal / (omega * mu * MU0);

      // tangential‐component
      cvec Et_contrib =  cross( (cvec)normal, (cvec)(beta * fval));
      cvec Ht_contrib = -cross( (cvec)normal, (cvec)(alpha * fval));

      // J,M currents have contributions from each RWG with opposite sign, thus -=
      Etot[triIndex][q] -= En_contrib + Et_contrib;
      Htot[triIndex][q] -= Hn_contrib + Ht_contrib;
    }
  }

  // Quadrature - find the scattering cross section
  double sigma = 0.0;
  for (size_t t = 0; t < Ntri; ++t) {
    rvec normal = vNormals[t];
    double area = vTri[t]->Area();
    
    // re-compute nodes for this triangle
    rvec r[3] = { vTri[t]->Node(0), vTri[t]->Node(1), vTri[t]->Node(2) };

    for (int q = 0; q < nQuad; ++q) {
      // map barycentrics -> point
      rvec s = xQuad[q][0] * r[0] + xQuad[q][1] * r[1] + xQuad[q][2] * r[2];

      cvec Einc = IncFieldE(s);
      cvec Hinc = IncFieldH(s);
      rvec Sinc = realVecHom( cross(Einc, conjVecHom(Hinc)) ) / 2.0;

      cvec E = Etot[t][q];
      cvec H = Htot[t][q];
      rvec S = realVecHom( cross(E, conjVecHom(H)) ) / 2.0;

      sigma -= wQuad[q] * area * dot(normal, S) / sqrt( dot( Sinc, Sinc ) );
    }
  }
  return sigma;
}

/**
 * @details Extinction cross section via interference of incident and scattered fields.
 */
double DomainHom3D::extinctionCS(blitz::Array<dcmplx,1>* solVector)
{
  // Gather unique triangles & normals
  std::vector<Triangle*> vTri;
  std::vector<rvec> vNormals;
  vTri.reserve(rwgFuns.size());
  vNormals.reserve(rwgFuns.size());

  for (RWGFun& rwg : rwgFuns) {
    Triangle* T = rwg.TrianglePtr();
    auto it = std::find(vTri.begin(), vTri.end(), T);
    if (it == vTri.end()) {
      vTri.push_back(T);
      vNormals.push_back(T->FBNormalVector(rwg.Side()));
    }
  }

  size_t Ntri = vTri.size();  // number of unique triangles

  // Quadrature rule on reference triangle
  int nQuad = dunavant::N5;    // number of quadrature points
  auto& xQuad = dunavant::x5;  // barycentrics [nQuad][3]
  auto& wQuad = dunavant::w5;  // weights[nQuad]

  // Allocate per-triangle, per-quad-point accumulators
  std::vector<std::vector<cvec>> Etot(Ntri, std::vector<cvec>(nQuad, cvec(0, 0, 0))),
                                 Htot(Ntri, std::vector<cvec>(nQuad, cvec(0, 0, 0)));

  int offset = solVector->shape()(0) / 2; // offset for β coefficients

  // Fill E, H normal and tangential components at each sample point
  for (RWGFun& rwg : rwgFuns) {
    // find which triangle index
    Triangle* T = rwg.TrianglePtr();
    size_t triIndex = std::distance(vTri.begin(),
                      std::find(vTri.begin(), vTri.end(), T));
    rvec normal = vNormals[triIndex];

    // RWG coefficients
    dcmplx alpha = (*solVector)(rwg.EdgeIndex());
    dcmplx beta  = (*solVector)(rwg.EdgeIndex() + offset);
    double divf  = rwg.RWGPreFactor();

    // triangle nodes
    rvec r[3] = { T->Node(0), T->Node(1), T->Node(2) };

    // loop quad points
    for (int q = 0; q < nQuad; ++q) {
      // physical sample point
      rvec s = xQuad[q][0] * r[0] + xQuad[q][1] * r[1] + xQuad[q][2] * r[2];

      // hRWG at s
      rvec fval = rwg.Evaluate(s);

      // normal‐component (constant over tri)
      cvec En_contrib = -I * alpha * divf * normal / (omega * epsilon * EPS0);
      cvec Hn_contrib = -I * beta * divf * normal / (omega * mu * MU0);

      // tangential‐component
      cvec Et_contrib =  cross( (cvec)normal, (cvec)(beta * fval));
      cvec Ht_contrib = -cross( (cvec)normal, (cvec)(alpha * fval));

      // J,M currents have contributions from each RWG with opposite sign, thus -=
      Etot[triIndex][q] -= En_contrib + Et_contrib;
      Htot[triIndex][q] -= Hn_contrib + Ht_contrib;
    }
  }

  // Quadrature - find the scattering cross section
  double sigma = 0.0;
  for (size_t t = 0; t < Ntri; ++t) {
    rvec normal = vNormals[t];
    double area = vTri[t]->Area();
    
    // re-compute nodes for this triangle
    rvec r[3] = { vTri[t]->Node(0), vTri[t]->Node(1), vTri[t]->Node(2) };

    for (int q = 0; q < nQuad; ++q) {
      // map barycentrics -> point
      rvec s = xQuad[q][0] * r[0] + xQuad[q][1] * r[1] + xQuad[q][2] * r[2];

      cvec Einc = IncFieldE(s);
      cvec Hinc = IncFieldH(s);
      rvec Sinc = realVecHom( cross(Einc, conjVecHom(Hinc)) ) / 2.0;

      cvec Escat = Etot[t][q] - Einc;
      cvec Hscat = Htot[t][q] - Hinc;
      rvec Sext  = realVecHom( cross(Einc, conjVecHom(Hscat)) + 
                               cross(Escat, conjVecHom(Hinc)) ) / 2.0;

      sigma -= wQuad[q] * area * dot(normal, Sext) / sqrt( dot( Sinc, Sinc ) );
    }
  }
  return sigma;
}

// -----------------------------
// Material properties & layered stubs
// -----------------------------
dcmplx DomainHom3D::Epsilon(rvec /*pos*/) { return epsilon; }

dcmplx DomainHom3D::Mu(rvec /*pos*/) { return mu; }

void DomainHom3D::newGrnFunLayered(std::vector<rvec>& /*posvec*/) {
  throw std::runtime_error("Not implemented for homogeneous domains.");
}

// -----------------------------
// Incident field factories
// -----------------------------
int DomainHom3D::AddPlaneWave(dcmplx wavelength, rvec propagationDirection,
                              cvec polarization) {
  IncidentField *ptr = new PlaneWave(wavelength, real(epsilon),
                                     propagationDirection, polarization);
  incidentFields.push_back(ptr);
  return 0;
}

int DomainHom3D::AddDipole(rvec location, cvec polarization) {
  IncidentField *ptr =
      new Dipole(location, polarization, vacuumWavelength, mu, grnFun);
  incidentFields.push_back(ptr);
  return 0;
}

int DomainHom3D::AddGaussian(double wavelength, double waist, rvec focal_point,
                             rvec propagationDirection, cvec polarization) {
  IncidentField *ptr =
      new Gaussian(wavelength, waist, focal_point, real(epsilon),
                   propagationDirection, polarization);
  incidentFields.push_back(ptr);
  return 0;
}
