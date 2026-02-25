/**
 * @file GreenF.cpp
 * @brief Implementation of the GreenF base class.
 */

#include "GreenF.h"

// -----------------------------
// Construction & lifetime
// -----------------------------
GreenF::GreenF(dcmplx wvl, dcmplx eps) : wavelength(wvl),
                                         k_B(2 * PI * csqrt(eps) / wvl) {}

GreenF::GreenF(dcmplx wvl, const std::vector<dcmplx>& eps)
    : wavelength(wvl), k_L(eps.size())  // Initialize k_L with the same size as eps
{
  for (size_t i = 0; i < eps.size(); ++i) {
    k_L[i] = 2.0 * PI * csqrt(eps[i]) / wvl;  // Compute k for each layer
  }
}

GreenF::GreenF(dcmplx wvl, dcmplx eps, dcmplx mu)
    : wavelength(wvl), k_B(2 * PI * csqrt(eps) * csqrt(mu) / wvl) {}

GreenF::GreenF(dcmplx wvl, 
               const std::vector<dcmplx>& eps, 
               const std::vector<dcmplx>& mu)
    : wavelength(wvl), k_L(eps.size())  // Initialize k_L with the same size as eps
{
  for (size_t i = 0; i < eps.size(); ++i) {
    k_L[i] = 2.0 * PI * csqrt(eps[i]) * csqrt(mu[i]) / wvl;  // Compute k_B for each layer
  }
}

GreenF::~GreenF() {}

// -----------------------------
// Accessors & utilities
// -----------------------------
dcmplx GreenF::WaveVec() { return k_B; }

/**
 * @brief Complex square root with controlled branch cut.
 * @details If the input is purely real, force the imaginary part to zero to be
 * explicit about the branch and then delegate to std::sqrt().
 */
dcmplx csqrt(dcmplx arg) {
  if (imag(arg) == 0) {
    dcmplx buf(real(arg), 0);
    arg = buf;
  }
  return sqrt(arg);
}

/**
 * @details Falls back to direct computation cerfc(z) when @c thereIsNoLookupTable is true
 * or when @p z lies outside the tabulated domain.
 */
dcmplx GreenF::EvErfc(dcmplx z) {
  // std::cout<<z<<"\n";
  if (thereIsNoLookupTable) {
    return cerfc(z);
  } else {
    double x((real(z) + t->maxRe) / t->incRe);
    double y((imag(z) + t->maxIm) / t->incIm);
    int i(x);
    int j(y);
    if ((x < 0.) || (y < 0.) || (i > t->sizeRe - 2) || (j > t->sizeIm - 2)) {
      return cerfc(z);
    } else {
      const double dx1 = x - i, dy1 = y - j;
      const double dx2 = 1. - dx1, dy2 = 1. - dy1;
      return (t->fun[i][j] * dx2 + t->fun[i + 1][j] * dx1) * dy2 +
             (t->fun[i][j + 1] * dx2 + t->fun[i + 1][j + 1] * dx1) * dy1;
    }
  }
}

int GreenF::InitLookupTable(LookupTableBin* tin) {
  t = tin;
  thereIsNoLookupTable = 0;
  return 0;
}

// -----------------------------
// Global switches
// -----------------------------
bool GreenF::NeedsAccurate = false;

void GreenF::EnableAccurate() {
  std::cout << "Using high-accuracy integration." << std::endl;
  NeedsAccurate = true;
}

bool GreenF::RequiresAccurate() { return NeedsAccurate; }

int GreenF::Etm = 1;

void GreenF::AssignEtm(int n) { Etm = n; }

// -----------------------------
// Layered wavenumber access
// -----------------------------
dcmplx GreenF::GetLayerWavenumber(size_t index) const {
  if (index >= k_L.size()) {
    throw std::out_of_range("Index out of bounds for k_L.");
  }
  return k_L[index];
}

std::vector<dcmplx> GreenF::GetLayersWavenumbers() const { return k_L; }
