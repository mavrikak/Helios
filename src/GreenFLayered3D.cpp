/**
 * @file GreenFLayered3D.cpp
 * @brief Implementation of the 3D Green's Function class in a layered background.
 */

#include "GreenFLayered3D.h"
#include <iostream>
#include <tuple>
#include "quadrature.h"
#include <algorithm>

/**
 * @brief Element-wise addition of two string–complex maps.
 * @param lhs Left-hand operand (copied before modification).
 * @param rhs Right-hand operand.
 * @return A map containing the sum of corresponding entries.
 * @details Missing keys in either map are treated as zero-valued.
 */
std::map<std::string, dcmplx>
operator+(std::map<std::string, dcmplx> lhs,
          std::map<std::string, dcmplx> const& rhs)
{
    // lhs now holds a copy of the first map's entries
    for (auto const& kv : rhs) {
        // kv.first is the key, kv.second is the value
        lhs[kv.first] += kv.second;  
        // if kv.first wasn't in lhs, lhs[kv.first] is default-constructed (0+0i)
    }
    return lhs;
}

/**
 * @brief Element-wise subtraction of two string–complex maps.
 * @param lhs Left-hand operand (copied before modification).
 * @param rhs Right-hand operand.
 * @return A map containing the difference of corresponding entries.
 * @details Missing keys in either map are treated as zero-valued.
 */
std::map<std::string, dcmplx>
operator-(std::map<std::string, dcmplx> lhs,
          std::map<std::string, dcmplx> const& rhs)
{
    // lhs now holds a copy of the first map's entries
    for (auto const& kv : rhs) {
        // kv.first is the key, kv.second is the value
        lhs[kv.first] -= kv.second;  
        // if kv.first wasn't in lhs, lhs[kv.first] is default-constructed (0+0i)
    }
    return lhs;
}

/**
 * @brief Scalar division of all values in a string–complex map.
 * @param lhs Input map (copied before division).
 * @param rhs Scalar divisor.
 * @return A map whose values are divided by @p rhs.
 * @details Performs complex-by-real division for each element.
 */
std::map<std::string,dcmplx>
operator/(std::map<std::string,dcmplx> lhs, double rhs)
{
  for(auto &kv : lhs)
    kv.second /= rhs;   // complex<double> /= double
  return lhs;
}

/**
 * @brief Creates a dyadic product of two 3D vectors.
 * @param a First vector.
 * @param b Second vector.
 * @return Dyadic product as a cdyad (3x3 matrix).
 */
cdyad dyadic(const cvec& a, const cvec& b) {
  cdyad result;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      result(i)(j) = a(i) * b(j);
  return result;
}

/**
 * @brief Computes the cross product of a vector with each row or column of a dyad.
 * @param v Vector to be crossed.
 * @param D Dyad (3x3 matrix).
 * @param rowWise If true, apply row-wise; if false, column-wise.
 * @return Resulting cdyad with cross products applied accordingly.
 */
cdyad vecCrossDyad(const cvec& v, const cdyad& D, bool rowWise = true){
  cdyad result;
  if (rowWise) {
    // Cross product with each row of the dyad
    for (int j = 0; j < 3; ++j)
      result(j) = cross(v, D(j));
  } else {
    // Cross product with each column of the dyad
    for (int j = 0; j < 3; ++j) {
      cvec col(D(0)[j], D(1)[j], D(2)[j]);
      cvec crossed = cross(v, col);
      result(0)[j] = crossed[0];
      result(1)[j] = crossed[1];
      result(2)[j] = crossed[2];
    }
  }

  return result;
}

extern "C" {
  /**
   * @brief External REGRIDPACK routine for evaluating a 2D interpolation function.
   *
   * Performs 2D interpolation of tabulated data @p f_colmajor over a grid
   * defined by @p x and @p y, at query points ( @p xq, @p yq ).
   *
   * @param mx Number of points along the x-direction of the input grid.
   * @param my Number of points along the y-direction of the input grid.
   * @param x Pointer to x-coordinates of the input grid (size mx).
   * @param y Pointer to y-coordinates of the input grid (size my).
   * @param f_colmajor Pointer to input function values in column-major order (size mx*my).
   * @param m Number of query points along x.
   * @param n Number of query points along y.
   * @param xq Pointer to query x-coordinates (size m).
   * @param yq Pointer to query y-coordinates (size n).
   * @param method Interpolation method selector (e.g., 0 = bilinear, 1 = bicubic).
   * @param fq_out Output array of interpolated values (size m*n).
   * @param ier Output error flag (0 if successful).
   */
  void rgrd2_eval_c(int mx, int my, const double* x, const double* y, 
                    const double* f_colmajor, int m, int n,
                    const double* xq, const double* yq, int method,
                    double* fq_out, int* ier);
}

/** 
 * @typedef intraTuple
 * @brief Tuple holding intra-layer coefficient data and metadata.
 */
using intraTuple = std::tuple<std::vector<dcmplx>, std::vector<int>, std::vector<int>, 
                              std::vector<double>, std::vector<double>>;

/**
 * @typedef interTuple
 * @brief Tuple holding inter-layer coefficient data and metadata.
 */
using interTuple = std::tuple<std::vector<double>, std::vector<double>, std::vector<double>, 
                              std::vector<dcmplx>>;

/** 
 * @brief Selector map from intra-layer interaction key to tuple accessor.
 *
 * Maps strings like "te", "tm", "tez", ... to a lambda returning the matching
 * @c intraTuple inside a @c CoeffsIntra object.
 */
std::map<std::string, std::function<intraTuple&(CoeffsIntra&)>> coeffSelectorIntra;

/**
 * @brief Selector map from inter-layer interaction key to tuple accessor.
 *
 * Maps strings like "te", "tm", "tez1", "tmz2", ... to a lambda returning the
 * matching @c interTuple inside a @c CoeffsInter object.
 */
std::map<std::string, std::function<interTuple&(CoeffsInter&)>> coeffSelectorInter;

/**
 * @brief Initialize the selector map for intra-layer coefficients.
 *
 * Populates @c coeffSelectorIntra with accessors keyed by terms:
 * "te", "tm", "tez", "tmz", "tezz", "tmzz", "tes", "tms", "ter", "tmr",
 * "terz", "tmrz". Each accessor returns a reference to the corresponding
 * @c intraTuple stored inside a @c CoeffsIntra instance.
 */
void intraMap() {
    coeffSelectorIntra["te"]   = std::function<intraTuple&(CoeffsIntra&)>([](CoeffsIntra& c) -> intraTuple& { return c.te; });
    coeffSelectorIntra["tm"]   = std::function<intraTuple&(CoeffsIntra&)>([](CoeffsIntra& c) -> intraTuple& { return c.tm; });
    coeffSelectorIntra["tez"]  = std::function<intraTuple&(CoeffsIntra&)>([](CoeffsIntra& c) -> intraTuple& { return c.tez; });
    coeffSelectorIntra["tmz"]  = std::function<intraTuple&(CoeffsIntra&)>([](CoeffsIntra& c) -> intraTuple& { return c.tmz; });
    coeffSelectorIntra["tezz"] = std::function<intraTuple&(CoeffsIntra&)>([](CoeffsIntra& c) -> intraTuple& { return c.tezz; });
    coeffSelectorIntra["tmzz"] = std::function<intraTuple&(CoeffsIntra&)>([](CoeffsIntra& c) -> intraTuple& { return c.tmzz; });
    coeffSelectorIntra["tes"]  = std::function<intraTuple&(CoeffsIntra&)>([](CoeffsIntra& c) -> intraTuple& { return c.tes; });
    coeffSelectorIntra["tms"]  = std::function<intraTuple&(CoeffsIntra&)>([](CoeffsIntra& c) -> intraTuple& { return c.tms; });
    coeffSelectorIntra["ter"]  = std::function<intraTuple&(CoeffsIntra&)>([](CoeffsIntra& c) -> intraTuple& { return c.ter; });
    coeffSelectorIntra["tmr"]  = std::function<intraTuple&(CoeffsIntra&)>([](CoeffsIntra& c) -> intraTuple& { return c.tmr; });
    coeffSelectorIntra["terz"] = std::function<intraTuple&(CoeffsIntra&)>([](CoeffsIntra& c) -> intraTuple& { return c.terz; });
    coeffSelectorIntra["tmrz"] = std::function<intraTuple&(CoeffsIntra&)>([](CoeffsIntra& c) -> intraTuple& { return c.tmrz; });
}

/**
 * @brief Initialize the selector map for inter-layer coefficients.
 *
 * Populates @c coeffSelectorInter with accessors keyed by terms:
 * "te", "tm", "tez1", "tez2", "tmz1", "tmz2", "tezz", "tmzz", "tes", "tms",
 * "ter", "tmr", "terz1", "terz2", "tmrz1", "tmrz2". Each accessor returns a
 * reference to the corresponding @c interTuple stored inside a @c CoeffsInter.
 */
void interMap() {
    coeffSelectorInter["te"]    = std::function<interTuple&(CoeffsInter&)>([](CoeffsInter& c) -> interTuple& { return c.te; });
    coeffSelectorInter["tm"]    = std::function<interTuple&(CoeffsInter&)>([](CoeffsInter& c) -> interTuple& { return c.tm; });
    coeffSelectorInter["tez1"]  = std::function<interTuple&(CoeffsInter&)>([](CoeffsInter& c) -> interTuple& { return c.tez1; });
    coeffSelectorInter["tez2"]  = std::function<interTuple&(CoeffsInter&)>([](CoeffsInter& c) -> interTuple& { return c.tez2; });
    coeffSelectorInter["tmz1"]  = std::function<interTuple&(CoeffsInter&)>([](CoeffsInter& c) -> interTuple& { return c.tmz1; });
    coeffSelectorInter["tmz2"]  = std::function<interTuple&(CoeffsInter&)>([](CoeffsInter& c) -> interTuple& { return c.tmz2; });
    coeffSelectorInter["tezz"]  = std::function<interTuple&(CoeffsInter&)>([](CoeffsInter& c) -> interTuple& { return c.tezz; });
    coeffSelectorInter["tmzz"]  = std::function<interTuple&(CoeffsInter&)>([](CoeffsInter& c) -> interTuple& { return c.tmzz; });
    coeffSelectorInter["tes"]   = std::function<interTuple&(CoeffsInter&)>([](CoeffsInter& c) -> interTuple& { return c.tes; });
    coeffSelectorInter["tms"]   = std::function<interTuple&(CoeffsInter&)>([](CoeffsInter& c) -> interTuple& { return c.tms; });
    coeffSelectorInter["ter"]   = std::function<interTuple&(CoeffsInter&)>([](CoeffsInter& c) -> interTuple& { return c.ter; });
    coeffSelectorInter["tmr"]   = std::function<interTuple&(CoeffsInter&)>([](CoeffsInter& c) -> interTuple& { return c.tmr; });
    coeffSelectorInter["terz1"] = std::function<interTuple&(CoeffsInter&)>([](CoeffsInter& c) -> interTuple& { return c.terz1; });
    coeffSelectorInter["terz2"] = std::function<interTuple&(CoeffsInter&)>([](CoeffsInter& c) -> interTuple& { return c.terz2; });
    coeffSelectorInter["tmrz1"] = std::function<interTuple&(CoeffsInter&)>([](CoeffsInter& c) -> interTuple& { return c.tmrz1; });
    coeffSelectorInter["tmrz2"] = std::function<interTuple&(CoeffsInter&)>([](CoeffsInter& c) -> interTuple& { return c.tmrz2; });
}

// Evaluates a single point using a 2D interpolation function.
double GreenFLayered3D::evaluateSinglePoint(
    double x, double y,
    std::tuple<std::vector<dcmplx>, std::vector<int>, std::vector<int>,
               std::vector<double>, std::vector<double>>& inCoeffs,
    std::string type)
{
  const std::vector<dcmplx>& F_flat_c = std::get<0>(inCoeffs);
  const std::vector<int>&     px      = std::get<1>(inCoeffs);
  const std::vector<int>&     py      = std::get<2>(inCoeffs);
  const std::vector<double>&  Xc      = std::get<3>(inCoeffs);
  const std::vector<double>&  Yc      = std::get<4>(inCoeffs);

  int mx = (type == "imag" ? py.size(), px.at(1) : px.at(0)); // both entries equal by construction
  int my = (type == "imag" ? py.at(1) : py.at(0));

  // Extract grids
  std::vector<double> x_grid(mx), y_grid(my);
  for (int i = 0; i < mx; ++i) x_grid[i] = Xc[i];
  for (int j = 0; j < my; ++j) y_grid[j] = Yc[j];

  // Clamp query to grid box
  auto clamp = [](double v, double lo, double hi) {
      return v < lo ? lo : (v > hi ? hi : v);
  };
  x = clamp(x, x_grid.front(), x_grid.back());
  y = clamp(y, y_grid.front(), y_grid.back());

  // Build the work array expected by rgrd2: F in column-major
  std::vector<double> F_flat(F_flat_c.size());
  if (type == "real") {
      std::transform(F_flat_c.begin(), F_flat_c.end(), F_flat.begin(),
                     [](const dcmplx& c){ return std::real(c); });
  } else {
      std::transform(F_flat_c.begin(), F_flat_c.end(), F_flat.begin(),
                     [](const dcmplx& c){ return std::imag(c); });
  }

  // Single-point interpolation
  const int m = 1, n = 1;
  const double xq = x, yq = y;
  double fq = 0.0;
  int ier = 0;

  // method: 1 => bicubic; change to 0 if you want strict bilinear.
  int method = 1;

  rgrd2_eval_c(mx, my, x_grid.data(), y_grid.data(), F_flat.data(),
               m, n, &xq, &yq, method, &fq, &ier);

  if (ier != 0) {
      std::cerr << "[regridpack] rgrd2_eval_c failed, ier=" << ier << "\n";
  }

  return fq;
}

// Evaluates a single point using a 3D interpolation function.
double GreenFLayered3D::evaluateSinglePoint(
    double x, double y, double z,
    std::tuple<std::vector<double>, std::vector<double>, std::vector<double>,
               std::vector<dcmplx>>& inCoeffs,
    std::string type)
{
    // 0 = linear in z (two 2D evals), 1 = Catmull–Rom cubic in z (four 2D evals)
    const int s_interpZMethod = 1;

    // Unpack tuple: Xv, Yv, Zv, Fc (complex volume, column-major: i + j*mx + k*mx*my)
    const std::vector<double>& Xv = std::get<0>(inCoeffs);
    const std::vector<double>& Yv = std::get<1>(inCoeffs);
    const std::vector<double>& Zv = std::get<2>(inCoeffs);
    const std::vector<dcmplx>& Fc = std::get<3>(inCoeffs);

    const int mx = static_cast<int>(Xv.size());
    const int my = static_cast<int>(Yv.size());
    const int mz = static_cast<int>(Zv.size());
    if (mx == 0 || my == 0 || mz == 0) return 0.0;

    // Build real-valued volume in the same layout
    std::vector<double> F(Fc.size());
    if (type == "real") {
        std::transform(Fc.begin(), Fc.end(), F.begin(),
                       [](const dcmplx& c){ return std::real(c); });
    } else { // assume "imag"
        std::transform(Fc.begin(), Fc.end(), F.begin(),
                       [](const dcmplx& c){ return std::imag(c); });
    }

    // Clamp (x,y,z) to grid box
    auto clamp = [](double v, double lo, double hi){ return v < lo ? lo : (v > hi ? hi : v); };
    const double xc = clamp(x, Xv.front(), Xv.back());
    const double yc = clamp(y, Yv.front(), Yv.back());
    const double zc = clamp(z, Zv.front(), Zv.back());

    // Find z-bracketing indices
    int k0 = 0, k1 = 0;
    if (mz == 1) {
        k0 = k1 = 0;
    } else {
        auto it = std::lower_bound(Zv.begin(), Zv.end(), zc);
        if (it == Zv.begin()) { k0 = k1 = 0; }
        else if (it == Zv.end()) { k0 = k1 = mz - 1; }
        else { k1 = static_cast<int>(it - Zv.begin()); k0 = k1 - 1; }
    }

    // Normalized parameter along z
    double t = 0.0;
    if (k1 != k0) {
        t = (zc - Zv[k0]) / (Zv[k1] - Zv[k0]);
        if (t < 0.0) t = 0.0;
        if (t > 1.0) t = 1.0;
    }

    // One-point 2D evaluation on a given z-slice
    const size_t slicePitch = static_cast<size_t>(mx) * static_cast<size_t>(my);
    const int m = 1, n = 1;
    const double xq = xc, yq = yc;
    const int method2d = 1; // 1 = bicubic, set 0 for bilinear

    auto evalSlice = [&](int kz) -> double {
        const double* Fk = F.data() + static_cast<size_t>(kz) * slicePitch;
        double fq = 0.0; int ier = 0;
        rgrd2_eval_c(mx, my, Xv.data(), Yv.data(), Fk, m, n, &xq, &yq, method2d, &fq, &ier);
        if (ier != 0) {
            std::cerr << "[regridpack] rgrd2_eval_c failed on z-slice " << kz
                      << " (ier=" << ier << ")\n";
        }
        return fq;
    };

    if (s_interpZMethod == 0 || mz == 1) {
        // ---- Linear in z (two slices) ----
        const double f0 = evalSlice(k0);
        const double f1 = evalSlice(k1);
        return (k0 == k1) ? f0 : ((1.0 - t) * f0 + t * f1);
    } else {
        // ---- Catmull–Rom cubic in z (four slices) ----
        auto catmullRom = [](double fm1, double f0, double f1, double f2, double tt) {
            const double t2 = tt * tt, t3 = t2 * tt;
            return 0.5 * ( (2.0*f0)
                         + (-fm1 + f1) * tt
                         + (2.0*fm1 - 5.0*f0 + 4.0*f1 - f2) * t2
                         + (-fm1 + 3.0*f0 - 3.0*f1 + f2) * t3 );
        };

        const int km1 = (k0 > 0)       ? (k0 - 1) : k0;
        const int kp1 = (k1 < mz - 1)  ? (k1 + 1) : k1;

        const double fm1 = evalSlice(km1);
        const double f0  = evalSlice(k0);
        const double f1  = evalSlice(k1);
        const double f2  = evalSlice(kp1);

        return (k0 == k1) ? f0 : catmullRom(fm1, f0, f1, f2, t);
    }
}

// Constructor with given medium properties.
GreenFLayered3D::GreenFLayered3D(dcmplx inWavelength,
                                 const std::vector<dcmplx>& inEpsilonMedium,
                                 const std::vector<double>& inZValsInterfaces,
                                 const std::vector<double>& inThickness,
                                 const std::vector<Grid>& inTabulationGrids,
                                 const bool enableLineInt, int inMidLayerIndex) 
    : GreenF(inWavelength, inEpsilonMedium),
      epsilon(inEpsilonMedium),
      zValsInterfaces(inZValsInterfaces),
      thickness(inThickness),
      tabulationGrids(inTabulationGrids),
      midLayerIndex(inMidLayerIndex) {
      for (size_t i = 0; i < epsilon.size(); i++)
      {
        mu.push_back(dcmplx(1.0, 0.0));
      }
      layeredUtils = new LayeredMediaUtils(k_L, epsilon, mu, zValsInterfaces, thickness, midLayerIndex);
      FillGreenFTable();
      enableLineIntegrals = enableLineInt;
      if (enableLineIntegrals) std::cout << "Line integrals activated." << std::endl;
      if (!enableLineIntegrals) std::cout << "Line integrals deactivated." << std::endl;
      }

// Constructor for magnetic media.
GreenFLayered3D::GreenFLayered3D(dcmplx inWavelength,
                                 const std::vector<dcmplx>& inEpsilonMedium,
                                 const std::vector<dcmplx>& inMuMedium,
                                 const std::vector<double>& inZValsInterfaces,
                                 const std::vector<double>& inThickness,
                                 const std::vector<Grid>& inTabulationGrids,
                                 const bool enableLineInt, int inMidLayerIndex) 
    : GreenF(inWavelength, inEpsilonMedium, inMuMedium),
      epsilon(inEpsilonMedium),
      mu(inMuMedium),
      zValsInterfaces(inZValsInterfaces),
      thickness(inThickness),
      tabulationGrids(inTabulationGrids),
      midLayerIndex(inMidLayerIndex) {
      layeredUtils = new LayeredMediaUtils(k_L, epsilon, mu, zValsInterfaces, thickness, midLayerIndex);
      intraMap();
      interMap();
      FillGreenFTable();
      enableLineIntegrals = enableLineInt;
      if (enableLineIntegrals) std::cout << "Line integrals activated." << std::endl;
      if (!enableLineIntegrals) std::cout << "Line integrals deactivated." << std::endl;
      }

/**
 * @details This function identifies the index of the layer in which a given point resides 
 * based on its z-coordinate relative to the predefined interface positions.
 */
int GreenFLayered3D::GetLayerIndex(rvec pos) const {
  double z = pos[2];  // Extract the z-component of the position
  int index = -1;
  const double eps = 1e-6; // Small tolerance for floating-point comparisons

  // Check if point is on an interface
  for (size_t i = 0; i < zValsInterfaces.size(); ++i) {
    if (std::abs(z - zValsInterfaces[i]) <= eps) {
      // Interface i separates layer i and layer i+1
      int l_lower = i;
      int l_upper = i + 1;
      index = std::abs(l_lower - midLayerIndex) < std::abs(l_upper - midLayerIndex) 
            ? l_lower : l_upper; break;
    }
  }

  for (size_t i = 0; i < zValsInterfaces.size() - 1; ++i) {
    if (z > zValsInterfaces[i] + eps && z < zValsInterfaces[i + 1] - eps) {
      index = i + 1; break;
    }
  }

  // If z is beyond the last interface, return the index of the last layer
  if (z > zValsInterfaces.back() + eps) {
    index = zValsInterfaces.size();
  }
  
  // If z is below the first interface, return the index of the first layer
  if (z < zValsInterfaces.front() - eps) {
    index = 0;
  }

  return index;
}

// Checks whether the given positions exceed a predefined threshold.
bool GreenFLayered3D::AboveThreshold(rvec r, rvec rp) {
  return norm(k_B) * dot(r - rp, r - rp) > 0.4;
}

// Checks whether the given triangles exceed a predefined threshold.
bool GreenFLayered3D::AboveThreshold(rvec r, Triangle* Tp) {
  // Barycenter of source triangle
  rvec rp = Tp->Center();

  // Layer indices of point and triangle
  int rIndex = GetLayerIndex(r), rpIndex = GetLayerIndex(rp);
  
  // Radial distance in xy‐plane
  double dx = r[0] - rp[0], dy = r[1] - rp[1], drho = std::sqrt(dx * dx + dy * dy);

  double dTot = 0.0; // Total distance = distance to interfaces + radial distance

  // Both points in the bottom layer (below first interface)
  if (rIndex == 0 && rpIndex == 0) {
    dTot = drho + std::fabs(r[2]  - zValsInterfaces[0]) +
                  std::fabs(rp[2] - zValsInterfaces[0]);
  }
  // Both points are in the top layer (above last interface)
  else if (rIndex == zValsInterfaces.size() && rpIndex == zValsInterfaces.size()) {
    dTot = drho + std::fabs(r[2]  - zValsInterfaces.back()) +
                  std::fabs(rp[2] - zValsInterfaces.back());
  }
  // Both points are in the same inner layer
  else if (rIndex == rpIndex && rIndex > 0 && rIndex < zValsInterfaces.size()) {
    // distance to lower interface
    double dBelow = std::fabs(r[2]  - zValsInterfaces[rIndex - 1]) +
                    std::fabs(rp[2] - zValsInterfaces[rIndex - 1]);
    // distance to upper interface
    double dAbove = std::fabs(r[2]  - zValsInterfaces[rIndex]) +
                    std::fabs(rp[2] - zValsInterfaces[rIndex]);
    dTot = drho + std::min(dBelow, dAbove);
  }
  // Point r is one layer above point rp
  else if (rIndex == rpIndex + 1) {
    dTot = drho + std::fabs(r[2]  - zValsInterfaces[rpIndex]) +
                  std::fabs(rp[2] - zValsInterfaces[rpIndex]);
  }
  // Point rp is one layer above point r
  else if (rpIndex == rIndex + 1) {
    dTot = drho + std::fabs(r[2]  - zValsInterfaces[rIndex]) +
                  std::fabs(rp[2] - zValsInterfaces[rIndex]);
  }
  // Other non-adjacent cases: no need for singularity treatment
  else {
    return true;
  }
  
  // Calculate the lengths of the edges of the source triangle
  double lengthsTp[3];
  for (int i = 0; i < 3; ++i) {
    lengthsTp[i] = sqrt( dot( Tp->Node(i % 3) - Tp->Node((i + 1) % 3),
                              Tp->Node(i % 3) - Tp->Node((i + 1) % 3) ) );
  }

  // Calculate the radius of the triangle's circumcircle
  double tauTp = Tp->Perimeter() / 2.;
  double rTp  = lengthsTp[0] * lengthsTp[1] * lengthsTp[2] / 
                ( 4. * sqrt( tauTp * (tauTp - lengthsTp[0]) *
                                     (tauTp - lengthsTp[1]) * 
                                     (tauTp - lengthsTp[2]) ) );
  
  // Check if the distance is greater than twice the circumradius
  return dTot > 2.0 * rTp;
}

// Checks whether the given triangles exceed a predefined threshold.
bool GreenFLayered3D::AboveThreshold(Triangle* T, Triangle* Tp) {
  // Barycenters of triangles
  rvec r = T->Center(), rp = Tp->Center();

  // Layer indices of barycenters/triangles
  int rIndex = GetLayerIndex(r), rpIndex = GetLayerIndex(rp);
  
  // Radial distance in xy‐plane
  double dx = r[0] - rp[0], dy = r[1] - rp[1], drho = std::sqrt(dx * dx + dy * dy);

  double dTot = 0.0; // Total distance = distance to interfaces + radial distance

  // Both points in the bottom layer (below first interface)
  if (rIndex == 0 && rpIndex == 0) {
    dTot = drho + std::fabs(r[2]  - zValsInterfaces[0]) +
                  std::fabs(rp[2] - zValsInterfaces[0]);
  }
  // Both points are in the top layer (above last interface)
  else if (rIndex == zValsInterfaces.size() && rpIndex == zValsInterfaces.size()) {
    dTot = drho + std::fabs(r[2]  - zValsInterfaces.back()) +
                  std::fabs(rp[2] - zValsInterfaces.back());
  }
  // Both points are in the same inner layer
  else if (rIndex == rpIndex && rIndex > 0 && rIndex < zValsInterfaces.size()) {
    // distance to lower interface
    double dBelow = std::fabs(r[2]  - zValsInterfaces[rIndex - 1]) +
                    std::fabs(rp[2] - zValsInterfaces[rIndex - 1]);
    // distance to upper interface
    double dAbove = std::fabs(r[2]  - zValsInterfaces[rIndex]) +
                    std::fabs(rp[2] - zValsInterfaces[rIndex]);
    dTot = drho + std::min(dBelow, dAbove);
  }
  // Point r is one layer above point rp
  else if (rIndex == rpIndex + 1) {
    dTot = drho + std::fabs(r[2]  - zValsInterfaces[rpIndex]) +
                  std::fabs(rp[2] - zValsInterfaces[rpIndex]);
  }
  // Point rp is one layer above point r
  else if (rpIndex == rIndex + 1) {
    dTot = drho + std::fabs(r[2]  - zValsInterfaces[rIndex]) +
                  std::fabs(rp[2] - zValsInterfaces[rIndex]);
  }
  // Other non-adjacent cases: no need for singularity treatment
  else {
    return true;
  }
  
  // Calculate the lengths of the edges of the triangles
  double lengthsT[3], lengthsTp[3];
  for (int i = 0; i < 3; ++i) {
    lengthsT[i]  = sqrt( dot( T->Node(i % 3) - T->Node((i + 1)  % 3),
                              T->Node(i % 3) - T->Node((i + 1)  % 3) ) );
    lengthsTp[i] = sqrt( dot( Tp->Node(i % 3) - Tp->Node((i + 1) % 3),
                              Tp->Node(i % 3) - Tp->Node((i + 1) % 3) ) );
  }

  // Calculate the radius of the triangles' circumcircles
  double tauT = T->Perimeter() / 2., tauTp = Tp->Perimeter() / 2.;
  double rT   = lengthsT[0] * lengthsT[1] * lengthsT[2] / 
                ( 4. * sqrt( tauT * (tauT - lengthsT[0]) *
                                    (tauT - lengthsT[1]) * 
                                    (tauT - lengthsT[2]) ) );
  double rTp  = lengthsTp[0] * lengthsTp[1] * lengthsTp[2] / 
                ( 4. * sqrt( tauTp * (tauTp - lengthsTp[0]) *
                                     (tauTp - lengthsTp[1]) * 
                                     (tauTp - lengthsTp[2]) ) );
  
  // Check if the distance is greater than twice the sum of the circumradii
  return dTot > 2.0 * (rT + rTp);
}

// Writes tabulation grids to a file.
void GreenFLayered3D::writeGridsToFile(const std::string& filename, const std::vector<Grid>& grids) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + filename);
    }

    for (size_t i = 0; i < grids.size(); ++i) {
        file << "Grid " << i + 1 << "\n";
        file << "(Source, Observation) Layers = (" << grids[i].i1 << "," << grids[i].i2 << ")" << "\n";
        file << "Radial grid: ";
        for (const auto& r : grids[i].r) {
            file << r << " ";
        }
        file << "\nZ1 grid: ";
        for (const auto& z1 : grids[i].z1) {
            file << z1 << " ";
        }
        if (!grids[i].z2.empty()) {
            file << "\nZ2 grid: ";
            for (const auto& z2 : grids[i].z2) {
                file << z2 << " ";
            }
        }
        file << "\n\n";
    }

    file.close();
}

// Evaluates the scalar Green's function in free space.
dcmplx GreenFLayered3D::Evaluate(rvec r, rvec rp) {
  rvec dif(r - rp);
  double dis(sqrt(dot(dif, dif)));
  return exp(I * k_B * dis) / (4. * PI * dis);
}

// Evaluates the dyadic Green's function in free space.
cdyad GreenFLayered3D::EvaluateDyadic(rvec r, rvec rp) {
  double dx(r(0) - rp(0));
  double dy(r(1) - rp(1));
  double dz(r(2) - rp(2));
  double rr(sqrt(dx * dx + dy * dy + dz * dz));
  dcmplx ff((I * k_B * rr - 1.) / (k_B * k_B * rr * rr));
  dcmplx nf((3. - 3. * I * k_B * rr - k_B * k_B * rr * rr) /
            (k_B * k_B * rr * rr * rr * rr));
  return cdyad(cvec(1. + ff + nf * dx * dx, nf * dx * dy, nf * dx * dz),
               cvec(nf * dy * dx, 1. + ff + nf * dy * dy, nf * dy * dz),
               cvec(nf * dz * dx, nf * dz * dy, 1. + ff + nf * dz * dz)) *
         exp(I * k_B * rr) / (4. * PI * rr);
}

// Evaluates the dyadic Green's function in a layered medium.
cdyad GreenFLayered3D::EvaluateDyadic(rvec r, rvec rp, std::string dyadicType)
{
  cdyad result( cvec( dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0) ),
                cvec( dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0) ),
                cvec( dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0) ) );
  
  /* Set up the interaction type (intraTB, intraIn, inter),
   * coefficient and layer indices (nIndex, mIndex).
   * -----------------------------------------------------------
   * intraTB: intra-layer interaction in top/bottom layer,
   * intraIn: intra-layer interaction in the rest of the layers,
   * inter: inter-layer interactions,
   * nIndex: observation layer, mIndex: source layer
   */
  std::string interactionType = "unknown";
  int coeffIndex = -1;
  int nIndex = GetLayerIndex(r);
  int mIndex = GetLayerIndex(rp);

  // Define the needed parameters for intra-layer and inter-layer interactions:
  // Source and observation layers are the same and are the top or bottom layer
  if ((nIndex == mIndex) && (nIndex == 0 || nIndex == epsilon.size() - 1)) {
    interactionType = "intraTB";
    auto it = std::find(intraTBIndex.begin(), intraTBIndex.end(), nIndex);
    if (it == intraTBIndex.end()) throw std::runtime_error("Invalid intraTB coeffIndex in EvaluateDyadic");
    coeffIndex = std::distance(intraTBIndex.begin(), it);
  }
  // Source and observation layers are the same but they are not the first or last layer
  else if ((nIndex == mIndex) && (nIndex != 0 && nIndex != epsilon.size() - 1)) {
    interactionType = "intraIn";
    auto it = std::find(intraPMIndex.begin(), intraPMIndex.end(), nIndex);
    if (it == intraPMIndex.end()) throw std::runtime_error("Invalid intraIn coeffIndex in EvaluateDyadic");
    coeffIndex = std::distance(intraPMIndex.begin(), it);
  }
  // Source and observation layers are different
  else if (nIndex != mIndex){
    interactionType = "inter";
    std::complex<int> target(nIndex, mIndex);
    auto it = std::find(interIndeces.begin(), interIndeces.end(), target);
    if (it == interIndeces.end()) throw std::runtime_error("Invalid inter coeffIndex in EvaluateDyadic");
    coeffIndex = std::distance(interIndeces.begin(), it);
  }
  else {
    std::cerr << "Something went wrong with the layers' indices during assembly evaluation" << std::endl;
  }

  std::vector<dcmplx> epsN(2), epsM(2), muN(2), muM(2);
  std::vector<dcmplx> kN(2), kM(2), kNM(2), kMN(2);
  std::vector<std::string> TE(2), TM(2);
  
  // Parameters for delta calculation
  epsN[0] = epsilon[nIndex];  epsM[0] = epsilon[mIndex];
  muN[0]  = mu[nIndex];       muM[0]  = mu[mIndex];
  kN[0]   = k_L[nIndex];      kM[0]   = k_L[mIndex];
  kNM[0]  = 2.0 * PI * csqrt(epsN[0]) * csqrt(muM[0]) / wavelength;
  kMN[0]  = 2.0 * PI * csqrt(epsM[0]) * csqrt(muN[0]) / wavelength;
  TE[0] = "te"; TM[0] = "tm";

  // Parameters for kappa calculation: duality principle (eps <-> mu, TE <-> TM)
  epsN[1] = mu[nIndex];       epsM[1] = mu[mIndex];
  muN[1]  = epsilon[nIndex];  muM[1]  = epsilon[mIndex];
  kN[1]   = k_L[nIndex];      kM[1]   = k_L[mIndex];
  kNM[1]  = 2.0 * PI * csqrt(epsN[1]) * csqrt(muM[1]) / wavelength;
  kMN[1]  = 2.0 * PI * csqrt(epsM[1]) * csqrt(muN[1]) / wavelength;
  TE[1] = "tm"; TM[1] = "te";
  
  // Calculate reflected Green's function quasistatic part and its derivatives
  double dsp(1e-4);                                                // Displacement
  rvec rXplusDX(r[0] + dsp, r[1] + 0.0, r[2] + 0.0);               // Observation point + dx
  rvec rYplusDY(r[0] + 0.0, r[1] + dsp, r[2] + 0.0);               // Observation point + dy
  rvec rZplusDZ(r[0] + 0.0, r[1] + 0.0, r[2] + dsp);               // Observation point + dz
  std::vector<rvec> obsPoints = {rXplusDX, rYplusDY, rZplusDZ};    // Displaced observation points
  rvec rpXplusDX(rp[0] + dsp, rp[1] + 0.0, rp[2] + 0.0);           // Source point + dx
  rvec rpYplusDY(rp[0] + 0.0, rp[1] + dsp, rp[2] + 0.0);           // Source point + dy
  rvec rpZplusDZ(rp[0] + 0.0, rp[1] + 0.0, rp[2] + dsp);           // Source point + dz
  std::vector<rvec> srcPoints = {rpXplusDX, rpYplusDY, rpZplusDZ}; // Displaced source points
  
  // dipole source position
  rvec spVec = rp;

  // Calculate the radial distance between source and observation points
  double rhoDist = sqrt((r[0] - spVec[0]) * (r[0] - spVec[0]) +
                        (r[1] - spVec[1]) * (r[1] - spVec[1]));
  if (rhoDist < 1e-10) rhoDist = 1e-10;

  // Calculate the unit vectors
  rvec unitRho = ( rvec( r[0] - spVec[0], r[1] - spVec[1], 0.0 ) ) / rhoDist;
  rvec unitChi(1.0, 0.0, 0.0);
  rvec unitPsi(0.0, 1.0, 0.0);
  rvec unitZeta(0.0, 0.0, 1.0);
  rvec unitPhi = cross(unitZeta, unitRho);

  // Convert source and observation coordinates to layered coordinates (for smooth evaluation)
  layeredCoords coordSmooth = 
        this->layeredUtils->cartesianToLayeredTab(r, spVec, interactionType);
  
  /* Evaluate smooth -- without quasistatic term -- tabulated Green's functions, 
   * namely: gTE, gTM, θz(gTE), θz(gTM), θzθz'(gTE), θzθz'(gTM), ...,
   * based on the interaction type (intra- or inter-layer) and their spatial derivatives.
   */
  std::map<std::string, dcmplx> gSmooth;
  std::vector<std::map<std::string, dcmplx>> gSmoothGradObs(obsPoints.size());
  std::vector<std::map<std::string, dcmplx>> gSmoothGradSrc(srcPoints.size());
  std::vector<std::vector<std::map<std::string, dcmplx>>> 
  gSmoothGradMixed(obsPoints.size(), 
                   std::vector<std::map<std::string, dcmplx>>(srcPoints.size()));
  gSmooth = (interactionType == "inter")
          ? this->evalGreenSmoothInter(coordSmooth.r, coordSmooth.Z1,
                                       coordSmooth.Z2, coeffIndex)
          : this->evalGreenSmoothIntra(interactionType, coordSmooth.r,
                                       coordSmooth.Z1, coordSmooth.Z2, coeffIndex);
  // derivatives calculation
  for (int i = 0; i < obsPoints.size(); ++i) {
    auto cSmoothObs = 
          this->layeredUtils->cartesianToLayeredTab(obsPoints[i], spVec, interactionType);
    auto cSmoothSrc = 
          this->layeredUtils->cartesianToLayeredTab(r, srcPoints[i], interactionType);
    auto gSmoothObs = (interactionType == "inter")
                    ? this->evalGreenSmoothInter(cSmoothObs.r, cSmoothObs.Z1,
                                                 cSmoothObs.Z2, coeffIndex)
                    : this->evalGreenSmoothIntra(interactionType, cSmoothObs.r,
                                                 cSmoothObs.Z1, cSmoothObs.Z2, coeffIndex);
    auto gSmoothSrc = (interactionType == "inter")
                    ? this->evalGreenSmoothInter(cSmoothSrc.r, cSmoothSrc.Z1,
                                                 cSmoothSrc.Z2, coeffIndex)
                    : this->evalGreenSmoothIntra(interactionType, cSmoothSrc.r,
                                                 cSmoothSrc.Z1, cSmoothSrc.Z2, coeffIndex);
    gSmoothGradObs[i] = (gSmoothObs - gSmooth) / dsp;
    gSmoothGradSrc[i] = (gSmoothSrc - gSmooth) / dsp;
    for (int j = 0; j < srcPoints.size(); ++j) {
      auto cSmoothSrcJ = 
            this->layeredUtils->cartesianToLayeredTab(r, srcPoints[j], interactionType);
      auto cSmoothMixed = 
            this->layeredUtils->cartesianToLayeredTab(obsPoints[i], srcPoints[j], 
                                                      interactionType);
      auto gSmoothSrcJ = (interactionType == "inter")
                       ? this->evalGreenSmoothInter(cSmoothSrcJ.r, cSmoothSrcJ.Z1,
                                                    cSmoothSrcJ.Z2, coeffIndex)
                       : this->evalGreenSmoothIntra(interactionType, cSmoothSrcJ.r,
                                                    cSmoothSrcJ.Z1, cSmoothSrcJ.Z2, 
                                                    coeffIndex);
      auto gSmoothMixed = (interactionType == "inter")
                        ? this->evalGreenSmoothInter(cSmoothMixed.r, cSmoothMixed.Z1,
                                                     cSmoothMixed.Z2, coeffIndex)
                        : this->evalGreenSmoothIntra(interactionType, cSmoothMixed.r,
                                                     cSmoothMixed.Z1, cSmoothMixed.Z2, 
                                                     coeffIndex);
      gSmoothGradMixed[i][j] = 
                      (gSmoothMixed - gSmoothObs - gSmoothSrcJ + gSmooth) / (dsp * dsp);
    }
  }

  /* Evaluate quasistatic terms for Green's functions, namely: gTE_qs, 
   * gTM_qs, θz(gTE_qs), θzz(gTM_qs), ..., based on the interaction type
   * (intra- or inter-layer) and their spatial derivatives.
   */
  std::map<std::string, dcmplx> gQS;
  std::vector<std::map<std::string, dcmplx>> gQSGradObs(obsPoints.size());
  std::vector<std::map<std::string, dcmplx>> gQSGradSrc(srcPoints.size());
  std::vector<std::vector<std::map<std::string, dcmplx>>> 
  gQSGradMixed(obsPoints.size(), 
               std::vector<std::map<std::string, dcmplx>>(srcPoints.size()));
  
  std::map<std::string, dcmplx> gTotal;
  std::vector<std::map<std::string, dcmplx>> gTotalGradObs(obsPoints.size());
  std::vector<std::map<std::string, dcmplx>> gTotalGradSrc(srcPoints.size());
  std::vector<std::vector<std::map<std::string, dcmplx>>> 
  gTotalGradMixed(obsPoints.size(), 
                  std::vector<std::map<std::string, dcmplx>>(srcPoints.size()));
  
  // Convert source and observation coordinates to layered coordinates (for quasistatic evaluation)
  auto coordQS = this->layeredUtils->cartesianToLayeredQS(r, spVec, interactionType);
  
  gQS = (interactionType == "inter")
      ? this->evalGreenQSInter(r, spVec, coordQS.r, coordQS.Z1, coordQS.R1)
      : this->evalGreenQSIntra(r, interactionType, coordQS.r, coordQS.Z1, 
                               coordQS.Z2, coordQS.R1, coordQS.R2);
  gTotal = gSmooth + gQS; /* Green's function = smooth + quasistatic (+ operator 
                           * is overloaded for this variable type) */
  
  // derivatives calculation
  for (int i = 0; i < obsPoints.size(); ++i) {
    auto cQSObs = 
          this->layeredUtils->cartesianToLayeredQS(obsPoints[i], spVec, interactionType);
    auto cQSSrc = 
          this->layeredUtils->cartesianToLayeredQS(r, srcPoints[i], interactionType);
    auto gQSObs = (interactionType == "inter")
                ? this->evalGreenQSInter(obsPoints[i], spVec, cQSObs.r, cQSObs.Z1, cQSObs.R1)
                : this->evalGreenQSIntra(obsPoints[i], interactionType, cQSObs.r, 
                                         cQSObs.Z1, cQSObs.Z2, cQSObs.R1, cQSObs.R2);
    auto gQSSrc = (interactionType == "inter")
                ? this->evalGreenQSInter(r, srcPoints[i], cQSSrc.r, cQSSrc.Z1, cQSSrc.R1)
                : this->evalGreenQSIntra(srcPoints[i], interactionType, cQSSrc.r, 
                                         cQSSrc.Z1, cQSSrc.Z2, cQSSrc.R1, cQSSrc.R2);
    gQSGradObs[i] = (gQSObs - gQS) / dsp;
    gTotalGradObs[i] = gSmoothGradObs[i] + gQSGradObs[i];
    gQSGradSrc[i] = (gQSSrc - gQS) / dsp;
    gTotalGradSrc[i] = gSmoothGradSrc[i] + gQSGradSrc[i];
    for (int j = 0; j < srcPoints.size(); ++j) {
      auto cQSSrcJ = 
          this->layeredUtils->cartesianToLayeredQS(r, srcPoints[j], interactionType);
      auto cQSMixed = 
          this->layeredUtils->cartesianToLayeredQS(obsPoints[i], srcPoints[j], 
                                                   interactionType);
      auto gQSSrcJ = (interactionType == "inter")
                   ? this->evalGreenQSInter(r, srcPoints[j], cQSSrcJ.r, 
                                            cQSSrcJ.Z1, cQSSrcJ.R1)
                   : this->evalGreenQSIntra(srcPoints[j], interactionType, cQSSrcJ.r, 
                                            cQSSrcJ.Z1, cQSSrcJ.Z2, cQSSrcJ.R1, cQSSrcJ.R2);
      auto gQSMixed = (interactionType == "inter")
                    ? this->evalGreenQSInter(obsPoints[i], srcPoints[j], cQSMixed.r, 
                                             cQSMixed.Z1, cQSMixed.R1)
                    : this->evalGreenQSIntra(obsPoints[i], interactionType, 
                                             cQSMixed.r, cQSMixed.Z1, cQSMixed.Z2, 
                                             cQSMixed.R1, cQSMixed.R2);
      gQSGradMixed[i][j] = (gQSMixed - gQSObs - gQSSrcJ + gQS) / (dsp * dsp);
      gTotalGradMixed[i][j] = gSmoothGradMixed[i][j] + gQSGradMixed[i][j];
    }
  }

  // build a cvec of the three partials
  auto gradGTotalObs = [&]( const std::string &key ) -> cvec {
    return cvec{
      gTotalGradObs[0].at(key),
      gTotalGradObs[1].at(key),
      gTotalGradObs[2].at(key)
    };
  };

  // build a cvec of the three partials
  auto gradGTotalSrc = [&]( const std::string &key ) -> cvec {
    return cvec{
      gTotalGradSrc[0].at(key),
      gTotalGradSrc[1].at(key),
      gTotalGradSrc[2].at(key)
    };
  };

  // build a cdyad of the nine mixed partials
  auto gradGTotalMixed = [&]( const std::string &key ) -> cdyad {
    return cdyad{
      { gTotalGradMixed[0][0].at(key), gTotalGradMixed[0][1].at(key), gTotalGradMixed[0][2].at(key) },
      { gTotalGradMixed[1][0].at(key), gTotalGradMixed[1][1].at(key), gTotalGradMixed[1][2].at(key) },
      { gTotalGradMixed[2][0].at(key), gTotalGradMixed[2][1].at(key), gTotalGradMixed[2][2].at(key) }
    };
  };

  // compute Delta or Kappa dyad
  if ( !(dyadicType == "normal" || dyadicType == "curl") ) {
    std::cout << std::endl << dyadicType << std::endl;
    throw std::runtime_error("Invalid dyadic type in GreenFLayered3D class (must be 'normal' or 'curl').");
  }
  if (dyadicType == "normal") {
    cdyad Delta = 
      ( gradGTotalMixed(TM[0] + "zz") / (kNM[0] * kNM[0]) - gradGTotalMixed(TE[0]) ) +
      ( gTotal[TM[0]] * kMN[0] * kMN[0] - gTotal[TE[0] + "zz"] ) * 
      dyadic( (cvec)(unitZeta), (cvec)(unitZeta) );
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        result(i)(j) += Delta(i)(j);
      }
    }
  }
  else if (dyadicType == "curl") {
    cdyad Kappa = 
      gTotal[TE[1] + "r"] * kM[1]  * kM[1]  * dyadic( (cvec)(unitPhi), (cvec)(unitZeta) ) -
      gTotal[TM[1] + "r"] * kMN[1] * kMN[1] * dyadic( (cvec)(unitZeta), (cvec)(unitPhi) );
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        result(i)(j) += Kappa(i)(j);
      }
    }
  }
  if (interactionType == "inter") {
    if(dyadicType == "normal") {
      cdyad Delta = dyadic( (cvec)(unitZeta),
                            ( gradGTotalSrc(TM[0] + "z2") * ( muN[0] / muM[0] ) + 
                              gradGTotalSrc(TE[0] + "z1") ) ) +
                    dyadic( ( gradGTotalObs(TM[0] + "z1") * (epsM[0] / epsN[0]) + 
                              gradGTotalObs(TE[0] + "z2") ), (cvec)(unitZeta) ) +
                    ( dyadic( (cvec)(unitRho), (cvec)(unitRho) ) + 
                      dyadic( (cvec)(unitPhi), (cvec)(unitPhi) ) ) * gTotal[TE[0] + "s"];
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          result(i)(j) += Delta(i)(j);
        }
      }
    }
    else if (dyadicType == "curl") {
      cdyad Kappa = vecCrossDyad( (cvec)(unitZeta), gradGTotalMixed(TM[1] + "z1") ) *
                    (epsM[1] / epsN[1]) + 
                    vecCrossDyad( (cvec)(unitZeta), gradGTotalMixed(TE[1] + "z2"), false );
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          result(i)(j) += Kappa(i)(j);
        }
      }
    }
  } 
  else {
    if(dyadicType == "normal") {
      cdyad Delta = dyadic( (cvec)(unitZeta),
                            ( gradGTotalSrc(TM[0] + "z") * ( muN[0] / muM[0] ) + 
                              gradGTotalSrc(TE[0] + "z") ) ) +
                    dyadic( ( gradGTotalObs(TM[0] + "z") * (epsM[0] / epsN[0]) + 
                              gradGTotalObs(TE[0] + "z") ), (cvec)(unitZeta) ) +
                    ( dyadic( (cvec)(unitChi), (cvec)(unitChi) ) + 
                      dyadic( (cvec)(unitPsi), (cvec)(unitPsi) ) ) * gTotal[TE[0] + "s"];
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          result(i)(j) += Delta(i)(j);
        }
      }
    }
    else if (dyadicType == "curl") {
      cdyad Kappa = vecCrossDyad( (cvec)(unitZeta), gradGTotalMixed(TM[1] + "z") ) *
                    (epsM[1] / epsN[1]) + 
                    vecCrossDyad( (cvec)(unitZeta), gradGTotalMixed(TE[1] + "z"), false );
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          result(i)(j) += Kappa(i)(j);
        }
      }
    }
  }
  return result;
}

// Smoothed Green's function to avoid singularity at r = rp
dcmplx GreenFLayered3D::Smoothed(rvec r, rvec rp) {
  rvec dif(r - rp);
  double dis(sqrt(dot(dif, dif)));
  dcmplx out(I * k_B / (4. * PI));
  dcmplx id(1., 0.);
  if (dis != 0.)
    out = 1. / (4. * PI) *
          ((exp(I * k_B * dis) - id) / dis + k_B * k_B * dis / 2.);
  return out;
}

// Half-smoothed Green's function to avoid singularity at r = rp
dcmplx GreenFLayered3D::halfSmoothed(rvec r, rvec rp) {
  rvec dif(r - rp);
  double dis(sqrt(dot(dif, dif)));
  dcmplx out(I * k_B / (4. * PI));
  dcmplx id(1., 0.);
  if (dis != 0.) out = 1. / (4. * PI) * ((exp(I * k_B * dis) - id) / dis);
  return out;
}

// Gradient of the Green's function in free space
cvec GreenFLayered3D::Gradient(rvec r, rvec rp) {
  rvec dif(r - rp);
  double dis(sqrt(dot(dif, dif)));
  cvec out(-exp(I * k_B * dis) / (4. * PI * dis * dis) * (-1. / dis + I * k_B) *
           dif);
  return out;
}

// Evaluate all: G, grad G, dyadic G
std::tuple<dcmplx, cvec, cdyad> GreenFLayered3D::EvaluateAll(rvec r, rvec rp) {
  rvec dif(r - rp);
  double dx(dif(0));
  double dy(dif(1));
  double dz(dif(2));
  double rr(sqrt(dx * dx + dy * dy + dz * dz));
  dcmplx t1(exp(I * k_B * rr) / (4. * PI * rr));
  dcmplx t2((I * k_B * rr - 1.) / (rr * rr));
  dcmplx t3((3. - 3. * I * k_B * rr - k_B * k_B * rr * rr) /
            (rr * rr * rr * rr));
  std::tuple<dcmplx, cvec, cdyad> ret;
  std::get<0>(ret) = t1;
  std::get<1>(ret) = -t1 * t2 * dif;
  std::get<2>(ret) =
      cdyad(cvec(t2 + t3 * dx * dx, t3 * dx * dy, t3 * dx * dz),
            cvec(t3 * dy * dx, t2 + t3 * dy * dy, t3 * dy * dz),
            cvec(t3 * dz * dx, t3 * dz * dy, t2 + t3 * dz * dz)) *
      t1;
  return ret;
}

// Gradient of the smoothed Green's function in free space
cvec GreenFLayered3D::GradientSmoothed(rvec r, rvec rp) {
  rvec dif(r - rp);
  double dis(sqrt(dot(dif, dif)));
  cvec out(0., 0., 0.);
  dcmplx id(1., 0.);
  if (dis != 0)
    out = -dif / (4. * PI * dis) *
          ((-exp(I * k_B * dis) + id) / (dis * dis) +
           I * k_B * exp(I * k_B * dis) / dis + k_B * k_B / 2.);
  return out;
}

// D term over RWG functions (not implemented)
dcmplx GreenFLayered3D::IntegrateD(RWGFun* /*f*/, RWGFun* /*fp*/) { exit(1); }

// K term over RWG functions (not implemented)
dcmplx GreenFLayered3D::IntegrateK(RWGFun* /*f*/, RWGFun* /*fp*/) { exit(1); }

// delta term: \f$\int_{Tp}\f$ G(r,·) f'(·) dS' evaluated at observation r.
cvec GreenFLayered3D::IntegrateDelta(rvec r, RWGFun* fp) 
{
  cvec out(0., 0., 0.);
  Triangle* Tp = fp->TrianglePtr();
  rvec rp[3];
  for (int i = 0; i < 3; i++) {
    rp[i] = Tp->Node(i);
  }
  double Ap = Tp->Area();
  rvec Tc(Tp->Center());
  bool needsSingCanc = !AboveThreshold(r, Tp);
  int nQuadrature;
  double(*xQuadrature)[3];
  double* wQuadrature;
  if (NeedsAccurate) {
    nQuadrature = needsSingCanc ? dunavant::N5 : dunavant::N1;
    xQuadrature = needsSingCanc ? dunavant::x5 : dunavant::x1;
    wQuadrature = needsSingCanc ? dunavant::w5 : dunavant::w1;
  } else {
    nQuadrature = dunavant::N1;
    xQuadrature = dunavant::x1;
    wQuadrature = dunavant::w1;
  }

  /* Set up the interaction type (intraTB, intraIn, inter),
   * coefficient and layer indices (nIndex, mIndex).
   * -----------------------------------------------------------
   * intraTB: intra-layer interaction in top/bottom layer,
   * intraIn: intra-layer interaction in the rest of the layers,
   * inter: inter-layer interactions,
   * nIndex: observation layer, mIndex: source layer
   */
  std::string interactionType = "unknown";
  int coeffIndex = -1;
  int nIndex = GetLayerIndex(r);
  int mIndex = GetLayerIndex(Tc);

  // Define the needed parameters for intra-layer and inter-layer interactions:
  // Source and observation layers are the same and are the top or bottom layer
  if ((nIndex == mIndex) && (nIndex == 0 || nIndex == epsilon.size() - 1)) {
    interactionType = "intraTB";
    auto it = std::find(intraTBIndex.begin(), intraTBIndex.end(), nIndex);
    if (it == intraTBIndex.end()) throw std::runtime_error("Invalid intraTB coeffIndex in IntegrateDelta");
    coeffIndex = std::distance(intraTBIndex.begin(), it);
  }
  // Source and observation layers are the same but they are not the first or last layer
  else if ((nIndex == mIndex) && (nIndex != 0 && nIndex != epsilon.size() - 1)) {
    interactionType = "intraIn";
    auto it = std::find(intraPMIndex.begin(), intraPMIndex.end(), nIndex);
    if (it == intraPMIndex.end()) throw std::runtime_error("Invalid intraIn coeffIndex in IntegrateDelta");
    coeffIndex = std::distance(intraPMIndex.begin(), it);
  }
  // Source and observation layers are different
  else if (nIndex != mIndex){
    interactionType = "inter";
    std::complex<int> target(mIndex, nIndex);
    auto it = std::find(interIndeces.begin(), interIndeces.end(), target);
    if (it == interIndeces.end()) throw std::runtime_error("Invalid inter coeffIndex in IntegrateDelta");
    coeffIndex = std::distance(interIndeces.begin(), it);
  }
  else {
    std::cerr << "Something went wrong with the layers' indices during assembly evaluation" << std::endl;
  }
  
  dcmplx epsN = epsilon[nIndex];  dcmplx epsM = epsilon[mIndex];
  dcmplx muN  = mu[nIndex];       dcmplx muM  = mu[mIndex];
  dcmplx kNM  = 2.0 * PI * csqrt(epsN) * csqrt(muM) / wavelength;
  dcmplx kMN  = 2.0 * PI * csqrt(epsM) * csqrt(muN) / wavelength;
  std::string TE = "te";          std::string TM = "tm";

  // Calculate reflected Green's function quasistatic part and its derivatives
  double dsp(1e-4);                                             // Displacement
  rvec rXplusDX(r[0] + dsp, r[1] + 0.0, r[2] + 0.0);            // Observation point + dx
  rvec rYplusDY(r[0] + 0.0, r[1] + dsp, r[2] + 0.0);            // Observation point + dy
  rvec rZplusDZ(r[0] + 0.0, r[1] + 0.0, r[2] + dsp);            // Observation point + dz
  std::vector<rvec> obsPoints = {rXplusDX, rYplusDY, rZplusDZ}; // Displaced observation points
  
  for (int i = 0; i < nQuadrature; ++i)
  {
    // quadrature point in the source triangle
    rvec spVec = rp[0] * xQuadrature[i][0] + rp[1] * xQuadrature[i][1] +
                 rp[2] * xQuadrature[i][2];

    double dAp = Ap * wQuadrature[i]; // weighted area of the observation triangle

    // Calculate the radial distance between source and observation points
    double rhoDist = sqrt((r[0] - spVec[0]) * (r[0] - spVec[0]) +
                          (r[1] - spVec[1]) * (r[1] - spVec[1]));
    if (rhoDist < 1e-10) rhoDist = 1e-10;

    // Calculate the unit vectors
    rvec unitRho = ( rvec( r[0] - spVec[0], r[1] - spVec[1], 0.0 ) ) / rhoDist;
    rvec unitZeta(0.0, 0.0, 1.0);
    rvec unitPhi = cross(unitZeta, unitRho);

    // Convert source and observation coordinates to layered coordinates (for smooth evaluation)
    layeredCoords coordSmooth = 
          this->layeredUtils->cartesianToLayeredTab(r, spVec, interactionType);
    
    /* Evaluate smooth -- without quasistatic term -- tabulated Green's functions, 
     * namely: gTE, gTM, θz(gTE), θz(gTM), θzθz'(gTE), θzθz'(gTM), ...,
     * based on the interaction type (intra- or inter-layer) and their spatial derivatives.
     */
    std::map<std::string, dcmplx> gSmooth;
    std::vector<std::map<std::string, dcmplx>> gSmoothGrad;
    gSmoothGrad.resize(obsPoints.size());
    gSmooth = (interactionType == "inter")
            ? this->evalGreenSmoothInter(coordSmooth.r, coordSmooth.Z1,
                                         coordSmooth.Z2, coeffIndex)
            : this->evalGreenSmoothIntra(interactionType, coordSmooth.r,
                                         coordSmooth.Z1, coordSmooth.Z2, coeffIndex);
    // derivatives calculation
    for (int i = 0; i < obsPoints.size(); ++i) {
      auto cSmoothtemp = this->layeredUtils->cartesianToLayeredTab(obsPoints[i], spVec,
                                                                   interactionType);
      auto gSmoothtemp = (interactionType == "inter")
                        ? this->evalGreenSmoothInter(cSmoothtemp.r, cSmoothtemp.Z1,
                                                     cSmoothtemp.Z2, coeffIndex)
                        : this->evalGreenSmoothIntra(interactionType, cSmoothtemp.r,
                                                     cSmoothtemp.Z1, cSmoothtemp.Z2, coeffIndex);
      gSmoothGrad[i] = (gSmoothtemp - gSmooth) / dsp;
    }

    /* Evaluate quasistatic terms for Green's functions, namely: gTE_qs, 
     * gTM_qs, θz(gTE_qs), θzz(gTM_qs), ..., based on the interaction type
     * (intra- or inter-layer) and their spatial derivatives.
     */
    std::map<std::string, dcmplx> gQS;                      // Quasistatic Green's functions
    std::vector<std::map<std::string, dcmplx>> gQSGrad;     // Quasistatic Green's functions gradient
    gQSGrad.resize(obsPoints.size());
    std::map<std::string, dcmplx> gTotal;                   // Total Green's functions
    std::vector<std::map<std::string, dcmplx>> gTotalGrad;  // Total Green's functions gradient
    gTotalGrad.resize(obsPoints.size());
    if (!needsSingCanc)
    {
      // Convert source and observation coordinates to layered coordinates (for quasistatic evaluation)
      auto coordQS = 
           this->layeredUtils->cartesianToLayeredQS(r, spVec, interactionType);
      
      // If no singularity cancellation is needed, use the typical quasistatic functions
      gQS = (interactionType == "inter")
          ? this->evalGreenQSInter(r, spVec, coordQS.r, coordQS.Z1, coordQS.R1)
          : this->evalGreenQSIntra(r, interactionType, coordQS.r, coordQS.Z1, 
                                   coordQS.Z2, coordQS.R1, coordQS.R2);
      
      gTotal = gSmooth + gQS; /* Green's function = smooth + quasistatic (+ operator 
                               * is overloaded for this variable type) */
      
      // derivatives calculation
      for (int i = 0; i < obsPoints.size(); ++i) {
        auto cQStemp = this->layeredUtils->cartesianToLayeredQS(obsPoints[i], spVec, interactionType);
        auto gQStemp = (interactionType=="inter")
                     ? evalGreenQSInter(obsPoints[i], spVec, cQStemp.r, cQStemp.Z1, cQStemp.R1)
                     : evalGreenQSIntra(obsPoints[i], interactionType, cQStemp.r, cQStemp.Z1,
                                        cQStemp.Z2, cQStemp.R1, cQStemp.R2);
        gQSGrad[i] = (gQStemp - gQS) / dsp;
        gTotalGrad[i] = gSmoothGrad[i] + gQSGrad[i];
      }
    } else {
      // If singularity cancellation is needed, use the analytical quasistatic functions
      gTotal = gSmooth; // Only smooth Green's functions are used
      for (int i = 0; i < obsPoints.size(); ++i) {
        gTotalGrad[i] = gSmoothGrad[i];
      }
    }

    // build a cvec of the three partials
    auto gradGTotal = [&]( const std::string &key ) -> cvec {
      return cvec{
        gTotalGrad[0].at(key),
        gTotalGrad[1].at(key),
        gTotalGrad[2].at(key)
      };
    };

    rvec fpTot    = fp->Evaluate(spVec);
    double fpRho  = dot(fpTot, unitRho);
    double fpPhi  = dot(fpTot, unitPhi);
    double fpZeta = dot(fpTot, unitZeta);
    double divfp  = fp->RWGPreFactor();

    // compute delta
    out += dAp * ( (-gradGTotal(TM + "zz") / (kNM * kNM) + gradGTotal(TE) ) * divfp  +
                   ( gTotal[TM] * kMN * kMN - gTotal[TE + "zz"]) * unitZeta * fpZeta +
                   (fpRho * unitRho + fpPhi * unitPhi) * gTotal[TE + "s"] );
    if (interactionType == "inter") {
      out -= 
      dAp * ( ( gTotal[TM + "z2"] * (muN / muM) + gTotal[TE + "z1"] ) * unitZeta * divfp -
              ( gradGTotal(TM + "z1") * (epsM / epsN) + gradGTotal(TE + "z2") ) * fpZeta );
    } else {
      out -= 
      dAp * ( ( gTotal[TM + "z"] * (muN / muM) + gTotal[TE + "z"] ) * unitZeta * divfp -
              ( gradGTotal(TM + "z") * (epsM / epsN) + gradGTotal(TE + "z") ) * fpZeta );
    }
  }
  if (needsSingCanc) {
    TriQuadPol triQuad(
      Tp,                 // pointer to source triangle
      r,                  // observation point
      gausslegendre::N3,  // number of Gauss-Legendre points for the 1D rule
      gausslegendre::x3,  // table of 1D Gauss-Legendre points
      gausslegendre::w3   // table of 1D Gauss-Legendre weights
    );
    // Get the quadrature points and weights
    std::vector<rvec> quadPoints = triQuad.getQuadPoints();
    std::vector<double> quadWeights = triQuad.getQuadWeights();
#ifdef DEBUG
    writeVectorToFile(quadPoints, quadWeights,
                     "quadPoints_rho_" + std::to_string(r[0]) + "_" + 
                                         std::to_string(r[1]) + "_" + 
                                         std::to_string(r[2]) + 
                      "_Tp_" + std::to_string((uintptr_t)Tp) + ".txt");
#endif
    // Loop over the quadrature points
    for (int idx = 0; idx < quadWeights.size(); ++idx)
    {
      // quadrature point in the source triangle
      rvec spVec = rp[0] * quadPoints[idx][0] + rp[1] * quadPoints[idx][1] + 
                   rp[2] * quadPoints[idx][2];
      double dW = quadWeights[idx];  // weights for the observation triangle
      if (dW < 1e-10) continue;
      // Calculate the radial distance between source and observation points
      double rhoDist = sqrt((r[0] - spVec[0]) * (r[0] - spVec[0]) +
                            (r[1] - spVec[1]) * (r[1] - spVec[1]));
      if (rhoDist < 1e-10) rhoDist = 1e-10;

      // Calculate the unit vectors
      rvec unitRho = ( rvec( r[0] - spVec[0], r[1] - spVec[1], 0.0 ) ) / rhoDist;
      rvec unitZeta(0.0, 0.0, 1.0);
      rvec unitPhi = cross(unitZeta, unitRho);

      /* Evaluate quasistatic terms for Green's functions, namely: gTE_qs, 
       * gTM_qs, θz(gTE_qs), θzz(gTM_qs), ..., based on the interaction type
       * (intra- or inter-layer) and their spatial derivatives.
       */
      std::map<std::string, dcmplx> gQS;                  // Quasistatic Green's functions
      std::vector<std::map<std::string, dcmplx>> gQSGrad; // Quasistatic Green's functions gradient
      gQSGrad.resize(obsPoints.size());

      // Convert source and observation coordinates to layered coordinates (for quasistatic evaluation)
      auto coordQS = layeredUtils->cartesianToLayeredQS(r, spVec, interactionType);

      // Since singularity cancellation is needed, use the modified quasistatic functions
      gQS = (interactionType=="inter")
          ? evalGreenQSInter(r, spVec, coordQS.r, coordQS.Z1, coordQS.R1, needsSingCanc)
          : evalGreenQSIntra(r, interactionType, coordQS.r, coordQS.Z1,
                             coordQS.Z2, coordQS.R1, coordQS.R2, needsSingCanc);
      
      // derivatives calculation
      for (int i = 0; i < obsPoints.size(); ++i) {
        auto cQStemp = layeredUtils->cartesianToLayeredQS(obsPoints[i], spVec, interactionType);
        auto gQStemp = (interactionType=="inter")
                     ? evalGreenQSInter(obsPoints[i], spVec, cQStemp.r, cQStemp.Z1,
                                        cQStemp.R1, needsSingCanc)
                     : evalGreenQSIntra(obsPoints[i], interactionType, cQStemp.r, cQStemp.Z1,
                                        cQStemp.Z2, cQStemp.R1, cQStemp.R2, needsSingCanc);
        gQSGrad[i] = (gQStemp - gQS) / dsp;
      }

      // build a cvec of the three partials
      auto gradGQS = [&]( const std::string &key ) -> cvec {
        return cvec{
          gQSGrad[0].at(key),
          gQSGrad[1].at(key),
          gQSGrad[2].at(key)
        };
      };

      rvec fpTot    = fp->Evaluate(spVec);
      double fpRho  = dot(fpTot, unitRho);
      double fpPhi  = dot(fpTot, unitPhi);
      double fpZeta = dot(fpTot, unitZeta);
      double divfp  = fp->RWGPreFactor();
      
      // compute DeltaKappa
      if (interactionType == "intraTB") {
        out += 
        dW * ( -( gQS[TM + "Refl"] + gQS[TE + "Refl"] ) * 
                ( gQS["z"] * unitZeta * divfp - gradGQS("z") * fpZeta ) -
                gQS[TM + "Refl"] * gradGQS("zz") * divfp / (kNM * kNM) -
                gQS[TE + "Refl"] * gQS["zz"] * unitZeta * fpZeta + 
                gQS[TE + "Refl"] * gQS["s"] * ( unitRho * fpRho + unitPhi * fpPhi ) );
      }
      else if (interactionType == "intraIn") {
        for (int k = 0; k < 2; ++k) {
          std::string intNum = std::to_string(k + 1); // lower (= 1) or upper (= 2) layer contribution
          out += 
          dW * ( -( gQS[TM + "Refl" + intNum] + gQS[TE + "Refl" + intNum] ) * 
                  ( gQS["z" + intNum] * unitZeta * divfp - gradGQS("z" + intNum) * fpZeta ) -
                  gQS[TM + "Refl" + intNum] * gradGQS("zz" + intNum) * divfp / (kNM * kNM) -
                  gQS[TE + "Refl" + intNum] * gQS["zz" + intNum] * unitZeta * fpZeta + 
                  gQS[TE + "Refl" + intNum] * 
                  gQS["s" + intNum] * ( unitRho * fpRho + unitPhi * fpPhi ) );
        }
      }
      else if (interactionType == "inter") {
        out +=
        dW * ( - gQS[TM + "Trans"] * gradGQS("zz") * divfp / (kNM * kNM) - 
                 gQS[TM + "Trans"] * gQS["z2"] * unitZeta * divfp * (muN / muM) - 
                 gQS[TE + "Trans"] * gQS["z1"] * unitZeta * divfp + 
                 gQS[TM + "Trans"] * gradGQS("z1") * fpZeta * (epsM / epsN) + 
                 gQS[TE + "Trans"] * gradGQS("z2") * fpZeta - 
                 gQS[TE + "Trans"] * gQS["zz"] * unitZeta * fpZeta + 
                 gQS[TE + "Trans"] * gQS["s"] * ( unitRho * fpRho + unitPhi * fpPhi ) );
      }
      else {
        throw std::runtime_error("Unknown interaction type: " + interactionType);
      }
    }
  }
  return out;
}

// kappa term: \f$\int_{Tp}\f$ (\f$\nabla \times \f$G(r,·))xf'(·) dS' evaluated at observation r.
cvec GreenFLayered3D::IntegrateKappa(rvec r, RWGFun* fp) {
  cvec out(0., 0., 0.);
  Triangle* Tp = fp->TrianglePtr();
  rvec rp[3];
  for (int i = 0; i < 3; i++) {
    rp[i] = Tp->Node(i);
  }
  double Ap = Tp->Area();
  rvec Tc(Tp->Center());
  bool needsSingCanc = !AboveThreshold(r, Tp);
  int nQuadrature;
  double(*xQuadrature)[3];
  double* wQuadrature;
  if (NeedsAccurate) {
    nQuadrature = needsSingCanc ? dunavant::N5 : dunavant::N1;
    xQuadrature = needsSingCanc ? dunavant::x5 : dunavant::x1;
    wQuadrature = needsSingCanc ? dunavant::w5 : dunavant::w1;
  } else {
    nQuadrature = dunavant::N1;
    xQuadrature = dunavant::x1;
    wQuadrature = dunavant::w1;
  }

  /* Set up the interaction type (intraTB, intraIn, inter),
   * coefficient and layer indices (nIndex, mIndex).
   * -----------------------------------------------------------
   * intraTB: intra-layer interaction in top/bottom layer,
   * intraIn: intra-layer interaction in the rest of the layers,
   * inter: inter-layer interactions,
   * nIndex: observation layer, mIndex: source layer
   */
  std::string interactionType = "unknown";
  int coeffIndex = -1;
  int nIndex = GetLayerIndex(r);
  int mIndex = GetLayerIndex(Tc);

  // Define the needed parameters for intra-layer and inter-layer interactions:
  // Source and observation layers are the same and are the top or bottom layer
  if ((nIndex == mIndex) && (nIndex == 0 || nIndex == epsilon.size() - 1)) {
    interactionType = "intraTB";
    auto it = std::find(intraTBIndex.begin(), intraTBIndex.end(), nIndex);
    if (it == intraTBIndex.end()) throw std::runtime_error("Invalid intraTB coeffIndex in IntegrateKappa");
    coeffIndex = std::distance(intraTBIndex.begin(), it);
  }
  // Source and observation layers are the same but they are not the first or last layer
  else if ((nIndex == mIndex) && (nIndex != 0 && nIndex != epsilon.size() - 1)) {
    interactionType = "intraIn";
    auto it = std::find(intraPMIndex.begin(), intraPMIndex.end(), nIndex);
    if (it == intraPMIndex.end()) throw std::runtime_error("Invalid intraIn coeffIndex in IntegrateKappa");
    coeffIndex = std::distance(intraPMIndex.begin(), it);
  }
  // Source and observation layers are different
  else if (nIndex != mIndex){
    interactionType = "inter";
    std::complex<int> target(mIndex, nIndex);
    auto it = std::find(interIndeces.begin(), interIndeces.end(), target);
    if (it == interIndeces.end()) throw std::runtime_error("Invalid inter coeffIndex in IntegrateKappa");
    coeffIndex = std::distance(interIndeces.begin(), it);
  }
  else {
    std::cerr << "Something went wrong with the layers' indices during assembly evaluation" << std::endl;
  }

  dcmplx epsN = mu[nIndex];       dcmplx epsM = mu[mIndex];
  dcmplx muN  = epsilon[nIndex];  dcmplx kM   = k_L[mIndex];
  dcmplx kMN  = 2.0 * PI * csqrt(epsM) * csqrt(muN) / wavelength;
  std::string TE = "tm";          std::string TM = "te";

  // Calculate reflected Green's function quasistatic part and its derivatives
  double dsp(1e-4);                                             // Displacement
  rvec rXplusDX(r[0] + dsp, r[1] + 0.0, r[2] + 0.0);            // Observation point + dx
  rvec rYplusDY(r[0] + 0.0, r[1] + dsp, r[2] + 0.0);            // Observation point + dy
  rvec rZplusDZ(r[0] + 0.0, r[1] + 0.0, r[2] + dsp);            // Observation point + dz
  std::vector<rvec> obsPoints = {rXplusDX, rYplusDY, rZplusDZ}; // Displaced observation points
  
  for (int i = 0; i < nQuadrature; ++i)
  {
    // quadrature point in the source triangle
    rvec spVec = rp[0] * xQuadrature[i][0] + rp[1] * xQuadrature[i][1] +
                 rp[2] * xQuadrature[i][2];

    double dAp = Ap * wQuadrature[i]; // weighted area of the observation triangle

    // Calculate the radial distance between source and observation points
    double rhoDist = sqrt((r[0] - spVec[0]) * (r[0] - spVec[0]) +
                          (r[1] - spVec[1]) * (r[1] - spVec[1]));
    if (rhoDist < 1e-10) rhoDist = 1e-10;

    // Calculate the unit vectors
    rvec unitRho = ( rvec( r[0] - spVec[0], r[1] - spVec[1], 0.0 ) ) / rhoDist;
    rvec unitZeta(0.0, 0.0, 1.0);
    rvec unitPhi = cross(unitZeta, unitRho);

    // Convert source and observation coordinates to layered coordinates (for smooth evaluation)
    layeredCoords coordSmooth = 
          this->layeredUtils->cartesianToLayeredTab(r, spVec, interactionType);
    
    /* Evaluate smooth -- without quasistatic term -- tabulated Green's functions, 
     * namely: gTE, gTM, θz(gTE), θz(gTM), θzθz'(gTE), θzθz'(gTM), ...,
     * based on the interaction type (intra- or inter-layer) and their spatial derivatives.
     */
    std::map<std::string, dcmplx> gSmooth;
    std::vector<std::map<std::string, dcmplx>> gSmoothGrad;
    gSmoothGrad.resize(obsPoints.size());
    gSmooth = (interactionType == "inter")
            ? this->evalGreenSmoothInter(coordSmooth.r, coordSmooth.Z1,
                                         coordSmooth.Z2, coeffIndex)
            : this->evalGreenSmoothIntra(interactionType, coordSmooth.r,
                                         coordSmooth.Z1, coordSmooth.Z2, coeffIndex);
    // derivatives calculation
    for (int i = 0; i < obsPoints.size(); ++i) {
      auto cSmoothtemp = this->layeredUtils->cartesianToLayeredTab(obsPoints[i], spVec,
                                                                   interactionType);
      auto gSmoothtemp = (interactionType == "inter")
                        ? this->evalGreenSmoothInter(cSmoothtemp.r, cSmoothtemp.Z1,
                                                     cSmoothtemp.Z2, coeffIndex)
                        : this->evalGreenSmoothIntra(interactionType, cSmoothtemp.r,
                                                     cSmoothtemp.Z1, cSmoothtemp.Z2, coeffIndex);
      gSmoothGrad[i] = (gSmoothtemp - gSmooth) / dsp;
    }

    /* Evaluate quasistatic terms for Green's functions, namely: gTE_qs, 
     * gTM_qs, θz(gTE_qs), θzz(gTM_qs), ..., based on the interaction type
     * (intra- or inter-layer) and their spatial derivatives.
     */
    std::map<std::string, dcmplx> gQS;                      // Quasistatic Green's functions
    std::vector<std::map<std::string, dcmplx>> gQSGrad;     // Quasistatic Green's functions gradient
    gQSGrad.resize(obsPoints.size());
    std::map<std::string, dcmplx> gTotal;                   // Total Green's functions
    std::vector<std::map<std::string, dcmplx>> gTotalGrad;  // Total Green's functions gradient
    gTotalGrad.resize(obsPoints.size());
    if (!needsSingCanc)
    {
      // Convert source and observation coordinates to layered coordinates (for quasistatic evaluation)
      auto coordQS = 
           this->layeredUtils->cartesianToLayeredQS(r, spVec, interactionType);
      
      // If no singularity cancellation is needed, use the typical quasistatic functions
      gQS = (interactionType == "inter")
          ? this->evalGreenQSInter(r, spVec, coordQS.r, coordQS.Z1, coordQS.R1)
          : this->evalGreenQSIntra(r, interactionType, coordQS.r, coordQS.Z1, 
                                   coordQS.Z2, coordQS.R1, coordQS.R2);
      
      gTotal = gSmooth + gQS; /* Green's function = smooth + quasistatic (+ operator 
                               * is overloaded for this variable type) */
      
      // derivatives calculation
      for (int i = 0; i < obsPoints.size(); ++i) {
        auto cQStemp = this->layeredUtils->cartesianToLayeredQS(obsPoints[i], spVec, interactionType);
        auto gQStemp = (interactionType=="inter")
                     ? evalGreenQSInter(obsPoints[i], spVec, cQStemp.r, cQStemp.Z1, cQStemp.R1)
                     : evalGreenQSIntra(obsPoints[i], interactionType, cQStemp.r, cQStemp.Z1,
                                        cQStemp.Z2, cQStemp.R1, cQStemp.R2);
        gQSGrad[i] = (gQStemp - gQS) / dsp;
        gTotalGrad[i] = gSmoothGrad[i] + gQSGrad[i];
      }
    } else {
      // If singularity cancellation is needed, use the analytical quasistatic functions
      gTotal = gSmooth; // Only smooth Green's functions are used
      for (int i = 0; i < obsPoints.size(); ++i) {
        gTotalGrad[i] = gSmoothGrad[i];
      }
    }

    // build a cvec of the three partials
    auto gradGTotal = [&]( const std::string &key ) -> cvec {
      return cvec{
        gTotalGrad[0].at(key),
        gTotalGrad[1].at(key),
        gTotalGrad[2].at(key)
      };
    };

    rvec fpTot    = fp->Evaluate(spVec);
    double fpRho  = dot(fpTot, unitRho);
    double fpPhi  = dot(fpTot, unitPhi);
    double fpZeta = dot(fpTot, unitZeta);
    double divfp  = fp->RWGPreFactor();

    // compute DeltaKappa
    out += dAp * ( gTotal[TE + "r"] * kM * kM * unitPhi *  fpZeta -
                   gTotal[TM + "r"] * kMN * kMN * unitZeta * fpPhi );
    if (interactionType == "inter") {
      out +=
      dAp * ( gTotal[TM + "rz1"] * (epsM / epsN) * unitPhi * fpRho / rhoDist -
              gradGTotal(TM + "rz1") * (epsM / epsN) * fpPhi -
              gTotal[TE + "rz2"] * unitPhi * divfp );
    } else {
      out +=
      dAp * ( gTotal[TM + "rz"] * (epsM / epsN) * unitPhi * fpRho / rhoDist -
              gradGTotal(TM + "rz") * (epsM / epsN) * fpPhi -
              gTotal[TE + "rz"] * unitPhi * divfp );
    }
  }
  if (needsSingCanc) {
    TriQuadPol triQuad(
      Tp,                 // pointer to source triangle
      r,                  // observation point
      gausslegendre::N3,  // number of Gauss-Legendre points for the 1D rule
      gausslegendre::x3,  // table of 1D Gauss-Legendre points
      gausslegendre::w3   // table of 1D Gauss-Legendre weights
    );
    // Get the quadrature points and weights
    std::vector<rvec> quadPoints = triQuad.getQuadPoints();
    std::vector<double> quadWeights = triQuad.getQuadWeights();
#ifdef DEBUG
    writeVectorToFile(quadPoints, quadWeights,
                     "quadPoints_rho_" + std::to_string(r[0]) + "_" + 
                                         std::to_string(r[1]) + "_" + 
                                         std::to_string(r[2]) + 
                      "_Tp_" + std::to_string((uintptr_t)Tp) + ".txt");
#endif
    // Loop over the quadrature points
    for (int idx = 0; idx < quadWeights.size(); ++idx)
    {
      // quadrature point in the source triangle
      rvec spVec = rp[0] * quadPoints[idx][0] + rp[1] * quadPoints[idx][1] + 
                   rp[2] * quadPoints[idx][2];
      double dW = quadWeights[idx];  // weights for the observation triangle
      if (dW < 1e-10) continue;
      // Calculate the radial distance between source and observation points
      double rhoDist = sqrt((r[0] - spVec[0]) * (r[0] - spVec[0]) +
                            (r[1] - spVec[1]) * (r[1] - spVec[1]));
      if (rhoDist < 1e-10) rhoDist = 1e-10;

      // Calculate the unit vectors
      rvec unitRho = ( rvec( r[0] - spVec[0], r[1] - spVec[1], 0.0 ) ) / rhoDist;
      rvec unitZeta(0.0, 0.0, 1.0);
      rvec unitPhi = cross(unitZeta, unitRho);

      /* Evaluate quasistatic terms for Green's functions, namely: gTE_qs, 
       * gTM_qs, θz(gTE_qs), θzz(gTM_qs), ..., based on the interaction type
       * (intra- or inter-layer) and their spatial derivatives.
       */
      std::map<std::string, dcmplx> gQS;                  // Quasistatic Green's functions
      std::vector<std::map<std::string, dcmplx>> gQSGrad; // Quasistatic Green's functions gradient
      gQSGrad.resize(obsPoints.size());

      // Convert source and observation coordinates to layered coordinates (for quasistatic evaluation)
      auto coordQS = layeredUtils->cartesianToLayeredQS(r, spVec, interactionType);

      // Since singularity cancellation is needed, use the modified quasistatic functions
      gQS = (interactionType=="inter")
          ? evalGreenQSInter(r, spVec, coordQS.r, coordQS.Z1, coordQS.R1, needsSingCanc)
          : evalGreenQSIntra(r, interactionType, coordQS.r, coordQS.Z1,
                             coordQS.Z2, coordQS.R1, coordQS.R2, needsSingCanc);
      
      // derivatives calculation
      for (int i = 0; i < obsPoints.size(); ++i) {
        auto cQStemp = layeredUtils->cartesianToLayeredQS(obsPoints[i], spVec, interactionType);
        auto gQStemp = (interactionType=="inter")
                     ? evalGreenQSInter(obsPoints[i], spVec, cQStemp.r, cQStemp.Z1,
                                        cQStemp.R1, needsSingCanc)
                     : evalGreenQSIntra(obsPoints[i], interactionType, cQStemp.r, cQStemp.Z1,
                                        cQStemp.Z2, cQStemp.R1, cQStemp.R2, needsSingCanc);
        gQSGrad[i] = (gQStemp - gQS) / dsp;
      }

      // build a cvec of the three partials
      auto gradGQS = [&]( const std::string &key ) -> cvec {
        return cvec{
          gQSGrad[0].at(key),
          gQSGrad[1].at(key),
          gQSGrad[2].at(key)
        };
      };

      rvec fpTot    = fp->Evaluate(spVec);
      double fpRho  = dot(fpTot, unitRho);
      double fpPhi  = dot(fpTot, unitPhi);
      double fpZeta = dot(fpTot, unitZeta);
      double divfp  = fp->RWGPreFactor();
      
      // compute DeltaKappa
      if (interactionType == "intraTB") {
        out += 
        dW * ( gQS[TM + "Refl"] * 
                ( gQS["rz"] * unitPhi * fpRho / rhoDist - gradGQS("rz") * fpPhi ) -
                gQS[TM + "Refl"] * gQS["r"] * unitZeta * fpPhi * kMN * kMN -
                gQS[TE + "Refl"] * gQS["rz"] * unitPhi * divfp + 
                gQS[TE + "Refl"] * gQS["r"] * unitPhi * fpZeta * kM * kM );
      }
      else if (interactionType == "intraIn") {
        for (int k = 0; k < 2; ++k) {
          std::string intNum = std::to_string(k + 1); // lower (= 1) or upper (= 2) layer contribution
          out += 
          dW * ( gQS[TM + "Refl" + intNum] * 
                  ( gQS["rz" + intNum] * unitPhi * fpRho / rhoDist - gradGQS("rz" + intNum) * fpPhi ) -
                  gQS[TM + "Refl" + intNum] * gQS["r" + intNum] * unitZeta * fpPhi * kMN * kMN -
                  gQS[TE + "Refl" + intNum] * gQS["rz" + intNum] * unitPhi * divfp + 
                  gQS[TE + "Refl" + intNum] * gQS["r" + intNum] * unitPhi * fpZeta * kM * kM );
        }
      }
      else if (interactionType == "inter") {
        out +=
        dW * ( gQS[TM + "Trans"] * 
                ( ( gQS["rz1"] * unitPhi * fpRho / rhoDist - gradGQS("rz1") * fpPhi ) * (epsM / epsN) - 
                  gQS["r"] * unitZeta * fpPhi * kMN * kMN ) - 
                gQS[TE + "Trans"] * 
                ( gQS["rz2"] * unitPhi * divfp - gQS["r"] * unitPhi * fpZeta * kM * kM) );
      }
      else {
        throw std::runtime_error("Unknown interaction type: " + interactionType);
      }
    }
  }
  return out;
}

std::pair<dcmplx, dcmplx> GreenFLayered3D::PairD(RWGFun* f, RWGFun* fp) {
  // Don't care about symmetry, this is more exact!
  dcmplx elem(IntegrateD(f, fp));
  if (f == fp) return std::pair<dcmplx, dcmplx>(elem, elem);
  dcmplx elem2(IntegrateD(fp, f));
  return std::pair<dcmplx, dcmplx>(elem, elem2);
}

std::pair<dcmplx, dcmplx> GreenFLayered3D::PairK(RWGFun* f, RWGFun* fp) {
  // Unfortunately, IntegrateK is not symmetric.
  dcmplx elem(IntegrateK(f, fp));
  if (f == fp) return std::pair<dcmplx, dcmplx>(elem, elem);
  dcmplx elem2(IntegrateK(fp, f));
  return std::pair<dcmplx, dcmplx>(elem, elem2);
}

// Triangle–triangle DK block assembly (same observation triangle group).
void GreenFLayered3D::SameTriDK(
    const std::vector<RWGFun*>& fvec, const std::vector<RWGFun*>& fpvec,
    std::vector<std::vector<std::pair<dcmplx, dcmplx>>>& DK)
{
  bool E_calc = true; /* Flag to check if DK_E or DK_H is calculated
                       * (according to DK input initialization)
                       */
  
  if (DK[0][0].first  == dcmplx(1.0, 0.0) && DK[0][0].second == dcmplx(0.0, 0.0)) {
    // Reset the entire matrix back to (0,0)
    for (auto &row : DK) {
      for (auto &entry : row) {
        entry.first  = dcmplx(0.0, 0.0);
        entry.second = dcmplx(0.0, 0.0);
      }
    }
    E_calc = false; // DK_E is not calculated, so we will calculate DK_H
  }
  else if (DK[0][0].first != dcmplx(0.0, 0.0) ||
           DK[0][0].second!= dcmplx(0.0, 0.0)) {
    throw std::runtime_error("DK matrix is not initialized correctly.");
  }

  Triangle *tPtr = fvec[0]->TrianglePtr(), *tpPtr = fpvec[0]->TrianglePtr();
  rvec r[3], rp[3];
  for (int i = 0; i < 3; i++) {
    r[i] = tPtr->Node(i);
    rp[i] = tpPtr->Node(i);
  }
  double A = tPtr->Area(), Ap = tpPtr->Area();
  rvec tCenter = tPtr->Center(), tpCenter = tpPtr->Center();
  int tristatus = tPtr->FindAdjacency(tpPtr);
  bool needsSingCanc = !AboveThreshold(tPtr, tpPtr) && tristatus > 0;
  int nQuadrature;
  double(*xQuadrature)[3];
  double* wQuadrature;
  if (NeedsAccurate) {
    nQuadrature = tristatus > 0 ? dunavant::N5 : dunavant::N1;
    xQuadrature = tristatus > 0 ? dunavant::x5 : dunavant::x1;
    wQuadrature = tristatus > 0 ? dunavant::w5 : dunavant::w1;
  } else {
    nQuadrature = dunavant::N1;
    xQuadrature = dunavant::x1;
    wQuadrature = dunavant::w1;
  }

  // Observation and source Dunavant quadrature points in each triangle
  std::vector<rvec> svec(nQuadrature), spvec(nQuadrature);
  for (int i = 0; i < nQuadrature; ++i) {
    svec[i] = r[0] * xQuadrature[i][0] + r[1] * xQuadrature[i][1] +
              r[2] * xQuadrature[i][2];
    spvec[i] = rp[0] * xQuadrature[i][0] + rp[1] * xQuadrature[i][1] +
               rp[2] * xQuadrature[i][2];
  }
  
  /* Set up the interaction type (intraTB, intraIn, inter),
   * coefficient and layer indices (nIndex, mIndex).
   * -----------------------------------------------------------
   * intraTB: intra-layer interaction in top/bottom layer,
   * intraIn: intra-layer interaction in the rest of the layers,
   * inter: inter-layer interactions,
   * nIndex: observation layer, mIndex: source layer
   */
  std::string interactionType = "unknown";
  int coeffIndex = -1;
  int nIndex = GetLayerIndex(tCenter);
  int mIndex = GetLayerIndex(tpCenter);

  // Define the needed parameters for intra-layer and inter-layer interactions:
  // Source and observation layers are the same and are the top or bottom layer
  if ((nIndex == mIndex) && (nIndex == 0 || nIndex == epsilon.size() - 1)) {
    interactionType = "intraTB";
    auto it = std::find(intraTBIndex.begin(), intraTBIndex.end(), nIndex);
    if (it == intraTBIndex.end()) throw std::runtime_error("Invalid intraTB coeffIndex in SameTriDK");
    coeffIndex = std::distance(intraTBIndex.begin(), it);
  }
  // Source and observation layers are the same but they are not the first or last layer
  else if ((nIndex == mIndex) && (nIndex != 0 && nIndex != epsilon.size() - 1)) {
    interactionType = "intraIn";
    auto it = std::find(intraPMIndex.begin(), intraPMIndex.end(), nIndex);
    if (it == intraPMIndex.end()) throw std::runtime_error("Invalid intraIn coeffIndex in SameTriDK");
    coeffIndex = std::distance(intraPMIndex.begin(), it);
  }
  // Source and observation layers are different
  else if (nIndex != mIndex){
    interactionType = "inter";
    std::complex<int> target(nIndex, mIndex);
    auto it = std::find(interIndeces.begin(), interIndeces.end(), target);
    if (it == interIndeces.end()) throw std::runtime_error("Invalid inter coeffIndex in SameTriDK");
    coeffIndex = std::distance(interIndeces.begin(), it);
  }
  else {
    std::cerr << "Something went wrong with the layers' indices during assembly evaluation" << std::endl;
  }

  // Check if line integrals on edges are needed
  edge e[3], ep[3];
  for (int i = 0; i < (int)fvec.size(); ++i)  { e[i]  = fvec[i]->Edge();  }
  for (int i = 0; i < (int)fpvec.size(); ++i) { ep[i] = fpvec[i]->Edge(); }

  bool needsTestLineInt[3]  = {false, false, false};  // Are test line integrals needed?
  bool needsBasisLineInt[3] = {false, false, false};  // Are basis line integrals needed?
  for (int j = 0; j < 3; j++) {
    // test/basis triangle edge j on interface i?
    rvec n1 = *(e[j].node1), np1 = *(ep[j].node1);
    rvec n2 = *(e[j].node2), np2 = *(ep[j].node2);
    for (int i = 0; i < (int)zValsInterfaces.size(); i++) {
      if (std::fabs(n1[2] - zValsInterfaces[i]) < 1e-10 &&
          std::fabs(n2[2] - zValsInterfaces[i]) < 1e-10 && enableLineIntegrals) {
        // activate test line integral on this edge iff the basis triangle's
        // layer is one of the two layers adjacent to interface i
        const bool testOn = (mIndex == i || mIndex == i + 1) && interIndeces.size() > 0;
        needsTestLineInt[j] = needsTestLineInt[j] || testOn;
        break;
      }
    }
    for (int i = 0; i < (int)zValsInterfaces.size(); i++) {
      if (std::fabs(np1[2] - zValsInterfaces[i]) < 1e-10 &&
          std::fabs(np2[2] - zValsInterfaces[i]) < 1e-10 && enableLineIntegrals) {
        // activate basis line integral on this edge iff the test triangle's
        // layer is one of the two layers adjacent to interface i
        const bool basisOn = (nIndex == i || nIndex == i + 1) && interIndeces.size() > 0;
        needsBasisLineInt[j] = needsBasisLineInt[j] || basisOn;
        break;
      }
    }
  }

#ifdef DEBUG
  std::vector<rvec> verts1 = {
    rvec(-1.02174878120422, 0.331926299259067, -24.9769061803818),
    rvec(-6.95718452334404, -4.40442040562630, -23.0),
    rvec(0.111776299308985, -8.00666809082031, -23.0)
  };
  std::vector<rvec> verts2 = {
    rvec(0.111776299308985, -8.00666809082031, -23.0),
    rvec(-6.95718452334404, -4.40442040562630, -23.0),
    rvec(-6.00274764001370, -11.7464363574982, -21.2364837527275)
  };

  // intra-layer interaction keys
  std::vector<std::string> intraKeys = {"te","tm","tez","tmz","tezz","tmzz",
                                        "tes","tms","ter","tmr","terz","tmrz"};
  // inter-layer interaction keys
  std::vector<std::string> interKeys = {"te","tm","tez1","tez2","tmz1","tmz2",
                                        "tezz","tmzz","tes","tms","ter","tmr",
                                        "terz1","terz2","tmrz1","tmrz2"};

  this->testSameTriDK(verts1, verts2, "te", "yes");
  for (int i = 1; i < interKeys.size(); ++i) {
    this->testSameTriDK(verts1, verts2, interKeys[i], "no");
  }
#endif

  dcmplx epsN(0.0, 0.0), epsM(0.0, 0.0), muN(0.0, 0.0), muM(0.0, 0.0);
  dcmplx kN(0.0, 0.0), kM(0.0, 0.0), kNM(0.0, 0.0), kMN(0.0, 0.0);
  std::string TE = "", TM = "";
  if (E_calc) { // Parameters for DK_E calculation
    epsN = epsilon[nIndex];  epsM = epsilon[mIndex];
    muN  = mu[nIndex];       muM  = mu[mIndex];
    kN   = k_L[nIndex];      kM   = k_L[mIndex];
    kNM  = 2.0 * PI * csqrt(epsN) * csqrt(muM) / wavelength;
    kMN  = 2.0 * PI * csqrt(epsM) * csqrt(muN) / wavelength;
    TE = "te"; TM = "tm";
  } else { // Parameters for DK_H calculation: duality principle (eps <-> mu, TE <-> TM)
    epsN = mu[nIndex];       epsM = mu[mIndex];
    muN  = epsilon[nIndex];  muM  = epsilon[mIndex];
    kN   = k_L[nIndex];      kM   = k_L[mIndex];
    kNM  = 2.0 * PI * csqrt(epsN) * csqrt(muM) / wavelength;
    kMN  = 2.0 * PI * csqrt(epsM) * csqrt(muN) / wavelength;
    TE = "tm"; TM = "te";
  }

  for (int i = 0; i < nQuadrature; ++i)
  {
    double dA = A * wQuadrature[i]; // Area of observation triangle
    for (int j = 0; j < nQuadrature; ++j)
    {
      double dAp = Ap * wQuadrature[j]; // Area of source triangle

      // Convert source and observation coordinates to layered coordinates (for smooth evaluation)
      layeredCoords coordSmooth = 
            this->layeredUtils->cartesianToLayeredTab(svec[i], spvec[j], interactionType);
      
      /* Evaluate smooth -- without quasistatic term -- tabulated Green's functions, 
       * namely: gTE, gTM, θz(gTE), θz(gTM), θzθz'(gTE), θzθz'(gTM), ...,
       * based on the interaction type (intra- or inter-layer).
       */
      auto gSmooth = (interactionType == "inter")
          ? this->evalGreenSmoothInter(coordSmooth.r, coordSmooth.Z1,
                                       coordSmooth.Z2, coeffIndex)
          : this->evalGreenSmoothIntra(interactionType, coordSmooth.r,
                                       coordSmooth.Z1, coordSmooth.Z2, coeffIndex);

      /* Evaluate quasistatic terms for Green's functions, namely: gTE_qs, 
       * gTM_qs, θz(gTE_qs), θzz(gTM_qs), ..., based on the interaction type
       * (intra- or inter-layer).
       */
      std::map<std::string, dcmplx> gQS;    // Quasistatic Green's functions
      std::map<std::string, dcmplx> gTotal; // Total Green's functions
      if (!needsSingCanc) {
        // Convert source and observation coordinates to layered coordinates (for quasistatic evaluation)
        layeredCoords coordQS = 
              this->layeredUtils->cartesianToLayeredQS(svec[i], spvec[j], interactionType);
        // If no singularity cancellation is needed, use the tabulated quasistatic functions
        gQS = (interactionType == "inter")
            ? this->evalGreenQSInter(svec[i], spvec[j], coordQS.r, coordQS.Z1, coordQS.R1)
            : this->evalGreenQSIntra(svec[i], interactionType, coordQS.r, coordQS.Z1, 
                                    coordQS.Z2, coordQS.R1, coordQS.R2);
        gTotal = gSmooth + gQS;  // Green's function = smooth + quasistatic (+ operator 
                                 // is overloaded for this variable type)
      } else {
        // If singularity cancellation is needed, use the analytical quasistatic functions
        gTotal = gSmooth; // Only smooth Green's functions are used
      }
      
      // Calculate the radial distance between source and observation points
      double rhoDist = sqrt((svec[i][0] - spvec[j][0]) * (svec[i][0] - spvec[j][0]) +
                            (svec[i][1] - spvec[j][1]) * (svec[i][1] - spvec[j][1]));
      
      if (rhoDist < 1e-10) rhoDist = 1e-10; // Avoid division by zero
      
      for (int k = 0; k < (int)fvec.size(); ++k) {
        RWGFun* f = fvec[k];
        rvec fs(f->Evaluate(svec[i]));
        double divf(f->RWGPreFactor());
        double fKterm = (fs[0] * (svec[i][1] - spvec[j][1]) - 
                         fs[1] * (svec[i][0] - spvec[j][0])) / rhoDist;
        for (int l = 0; l < (int)fpvec.size(); ++l) {
          RWGFun* fp = fpvec[l];
          rvec fpsp(fp->Evaluate(spvec[j]));
          double divfp(fp->RWGPreFactor());
          double fpKterm = (fpsp[0] * (svec[i][1] - spvec[j][1]) - 
                            fpsp[1] * (svec[i][0] - spvec[j][0])) / rhoDist;
          DK[k][l].first += 
          dA * dAp * (divf * (gTotal[TM + "zz"] / (kNM * kNM) - gTotal[TE]) * divfp  +
                      fs[2] * (gTotal[TM] * kMN * kMN - gTotal[TE + "zz"]) * fpsp[2] +
                      fs[0] * gTotal[TE + "s"] * fpsp[0] + fs[1] * gTotal[TE + "s"] * fpsp[1]);
          DK[k][l].second -=  // Opposite sign compared to Chew's formula
          dA * dAp * (fKterm * gTotal[TE + "r"] * kM * kM * fpsp[2] -
                      fs[2] * gTotal[TM + "r"] * kMN * kMN * fpKterm);
          if (interactionType == "inter"){
            DK[k][l].first -= 
            dA * dAp * (fs[2] * (gTotal[TM + "z2"] * (muN / muM) + gTotal[TE + "z1"]) * divfp +
                        divf * (gTotal[TM + "z1"] * (epsM / epsN) + gTotal[TE + "z2"]) * fpsp[2]);
            DK[k][l].second -=  // Opposite sign compared to Chew's formula 
            dA * dAp * (divf * gTotal[TM + "rz1"] * (epsM / epsN) * fpKterm -
                        fKterm * gTotal[TE + "rz2"] * divfp);
          } else {
            DK[k][l].first -= 
            dA * dAp * (fs[2] * (gTotal[TM + "z"] * (muN / muM) + gTotal[TE + "z"]) * divfp +
                        divf * (gTotal[TM + "z"] * (epsM / epsN) + gTotal[TE + "z"]) * fpsp[2]);
            DK[k][l].second -=  // Opposite sign compared to Chew's formula 
            dA * dAp * (divf * gTotal[TM + "rz"] * (epsM / epsN) * fpKterm -
                        fKterm * gTotal[TE + "rz"] * divfp);
          }
        }
      }
    }
  }
  if (needsSingCanc) {
    // 1) build the Taylor-Duffy quadrature on this triangle pair
    SingCancellation singCanc(
      tristatus,         // 1=vertex, 2=edge, 3=facet
      tPtr, tpPtr,       // pointers to the triangles
      gausslegendre::N3, // number of Gauss-Legendre points for the 1D rule
      gausslegendre::x3, // table of 1D Gauss-Legendre points
      gausslegendre::w3  // table of 1D Gauss-Legendre weights
    );
    // grab the two maps of barycentric coords plus the combined weights
    std::vector<std::vector<double>> xi = singCanc.getPointsMap(0);  // 9 x 9 x numOfIntegrals x 3
    std::vector<std::vector<double>> eta = singCanc.getPointsMap(1); // 9 x 9 x numOfIntegrals x 3
    std::vector<double> weight = singCanc.getWeights();              // 9 x 9 x numOfIntegrals

    // 2) loop over every quadrature point
    int Nsc = weight.size();
    for (int idx = 0; idx < Nsc; ++idx) {
      // reconstruct observation/source points in 3D
      rvec sVec  = r[0] * xi[idx][0] + r[1] * xi[idx][1] + r[2] * xi[idx][2];
      rvec spVec = rp[0] * eta[idx][0] + rp[1] * eta[idx][1] + rp[2] * eta[idx][2];
      // the weight on this subpatch already includes jacobians and weights
      double dW = weight[idx] * A * Ap;

      // quasi‐static Green's functions
      auto coordQS = layeredUtils->cartesianToLayeredQS(sVec, spVec, interactionType);
      auto gQS = (interactionType=="inter")
               ? evalGreenQSInter(sVec, spVec, coordQS.r, coordQS.Z1, coordQS.R1)
               : evalGreenQSIntra(sVec, interactionType, coordQS.r, coordQS.Z1,
                                  coordQS.Z2, coordQS.R1, coordQS.R2);

      // Calculate the radial distance between source and observation points
      double rhoDist = sqrt((sVec[0] - spVec[0]) * (sVec[0] - spVec[0]) +
                            (sVec[1] - spVec[1]) * (sVec[1] - spVec[1]));
      if (rhoDist < 1e-10) rhoDist = 1e-10;

      for (int k = 0; k < (int)fvec.size(); ++k) {
        RWGFun* f = fvec[k];
        rvec  fs = f->Evaluate(sVec);
        double divf = f->RWGPreFactor();
        double fKterm = (fs[0] * (sVec[1] - spVec[1]) - 
                         fs[1] * (sVec[0] - spVec[0])) / rhoDist;
        for (int l = 0; l < (int)fpvec.size(); ++l) {
          RWGFun* fp = fpvec[l];
          rvec  fpsp = fp->Evaluate(spVec);
          double divfp = fp->RWGPreFactor();
          double fpKterm = (fpsp[0] * (sVec[1] - spVec[1]) - 
                            fpsp[1] * (sVec[0] - spVec[0])) / rhoDist;

          // Assemble exactly as in the non-singular smooth case, but weight by dW:
          DK[k][l].first += 
          dW * (divf * (gQS[TM + "zz"] / (kNM * kNM) - gQS[TE]) * divfp  +
                fs[2] * (gQS[TM] * kMN * kMN - gQS[TE + "zz"]) * fpsp[2] +
                fs[0] * gQS[TE + "s"] * fpsp[0] + fs[1] * gQS[TE + "s"] * fpsp[1]);
          DK[k][l].second -=  // Opposite sign compared to Chew's formula
          dW * (fKterm * gQS[TE + "r"] * kM * kM * fpsp[2] -
                fs[2] * gQS[TM + "r"] * kMN * kMN * fpKterm);
          if (interactionType == "inter"){
            DK[k][l].first -= 
            dW * (fs[2] * (gQS[TM + "z2"] * (muN / muM) + gQS[TE + "z1"]) * divfp +
                  divf * (gQS[TM + "z1"] * (epsM / epsN) + gQS[TE + "z2"]) * fpsp[2]);
            DK[k][l].second -=  // Opposite sign compared to Chew's formula 
            dW * (divf * gQS[TM + "rz1"] * (epsM / epsN) * fpKterm -
                  fKterm * gQS[TE + "rz2"] * divfp);
          } else {
            DK[k][l].first -= 
            dW * (fs[2] * (gQS[TM + "z"] * (muN / muM) + gQS[TE + "z"]) * divfp +
                  divf * (gQS[TM + "z"] * (epsM / epsN) + gQS[TE + "z"]) * fpsp[2]);
            DK[k][l].second -=  // Opposite sign compared to Chew's formula 
            dW * (divf * gQS[TM + "rz"] * (epsM / epsN) * fpKterm -
                  fKterm * gQS[TE + "rz"] * divfp);
          }
        }
      }
    }
  }
  if ( (needsTestLineInt[0] || needsTestLineInt[1] || needsTestLineInt[2]) &&
      !(needsTestLineInt[0] && needsTestLineInt[1] && needsTestLineInt[2]) ) {
    // Line integral for test functions
    for (int i = 0; i < 3; i++) {
      if (!needsTestLineInt[i]) continue;
      rvec n1 = *(e[i].node1), n2 = *(e[i].node2);
      rvec n = cross((rvec)(n1 - tCenter), (rvec)(n2 - tCenter));
      double L = sqrt(dot(n2 - n1, n2 - n1));

      // Outward normal on edge;
      rvec m = cross((rvec)(n2 - n1), n);
      m = m / sqrt(dot(m, m));

      // Line integral
      for (int j = 0; j < gausslegendre::N4; j++) {
        rvec s = (n1 - 1e-6 * m) * gausslegendre::x4[j][0] + 
                 (n2 - 1e-6 * m) * gausslegendre::x4[j][1];
        double dL = L * gausslegendre::w4[j];
        for (int q = 0; q < nQuadrature; ++q) {
          double dAp = Ap * wQuadrature[q];
          // Convert source and observation coordinates to layered coordinates (for smooth evaluation)
          layeredCoords coordSmooth = 
                this->layeredUtils->cartesianToLayeredTab(s, spvec[q], interactionType);
          
          /* Evaluate smooth -- without quasistatic term -- tabulated Green's functions, 
          * namely: gTE, gTM, θz(gTE), θz(gTM), θzθz'(gTE), θzθz'(gTM), ...,
          * based on the interaction type (intra- or inter-layer).
          */
          auto gSmooth = (interactionType == "inter")
              ? this->evalGreenSmoothInter(coordSmooth.r, coordSmooth.Z1,
                                           coordSmooth.Z2, coeffIndex)
              : this->evalGreenSmoothIntra(interactionType, coordSmooth.r,
                                           coordSmooth.Z1, coordSmooth.Z2, coeffIndex);

          /* Evaluate quasistatic terms for Green's functions, namely: gTE_qs, 
          * gTM_qs, θz(gTE_qs), θzz(gTM_qs), ..., based on the interaction type
          * (intra- or inter-layer).
          */
          std::map<std::string, dcmplx> gQS;    // Quasistatic Green's functions
          std::map<std::string, dcmplx> gTotal; // Total Green's functions
          
          // Convert source and observation coordinates to layered coordinates (for quasistatic evaluation)
          layeredCoords coordQS = 
                this->layeredUtils->cartesianToLayeredQS(s, spvec[q], interactionType);
          // If no singularity cancellation is needed, use the tabulated quasistatic functions
          gQS = (interactionType == "inter")
              ? this->evalGreenQSInter(s, spvec[q], coordQS.r, coordQS.Z1, coordQS.R1)
              : this->evalGreenQSIntra(s, interactionType, coordQS.r, coordQS.Z1, 
                                       coordQS.Z2, coordQS.R1, coordQS.R2);
          gTotal = gSmooth + gQS; // Green's function = smooth + quasistatic (+ operator 
                                  // is overloaded for this variable type)
          
          // Calculate the radial distance between source and observation points
          double rhoDist = sqrt((s[0] - spvec[q][0]) * (s[0] - spvec[q][0]) +
                                (s[1] - spvec[q][1]) * (s[1] - spvec[q][1]));
          
          if (rhoDist < 1e-10) rhoDist = 1e-10; // Avoid division by zero
          for (int k = 0; k < (int)fvec.size(); ++k) {
            RWGFun* f = fvec[k];
            rvec fs(f->Evaluate(s));
            for (int l = 0; l < (int)fpvec.size(); ++l) {
              RWGFun* fp = fpvec[l];
              rvec fpsp(fp->Evaluate(spvec[q]));
              double fpKterm = (fpsp[0] * (s[1] - spvec[q][1]) - 
                                fpsp[1] * (s[0] - spvec[q][0])) / rhoDist;
              if (interactionType == "inter"){
                DK[k][l].second +=  // Opposite sign compared to Chew's formula 
                (epsM / epsN) * dL * dAp * dot(fs, m) * fpKterm * gTotal[TM + "rz1"];
              } else {
                DK[k][l].second +=  // Opposite sign compared to Chew's formula 
                (epsM / epsN) * dL * dAp * dot(fs, m) * fpKterm * gTotal[TM + "rz"];
              }
            }
          }
        }
      }
    }
  }
  if ( (needsBasisLineInt[0] || needsBasisLineInt[1] || needsBasisLineInt[2]) &&
      !(needsBasisLineInt[0] && needsBasisLineInt[1] && needsBasisLineInt[2]) ) {
    // Line integral for basis functions
    for (int i = 0; i < 3; i++) {
      if (!needsBasisLineInt[i]) continue;
      rvec np1 = *(ep[i].node1), np2 = *(ep[i].node2);
      rvec np = cross((rvec)(np1 - tpCenter), (rvec)(np2 - tpCenter));
      double Lp = sqrt(dot(np2 - np1, np2 - np1));

      // Outward normal on edge;
      rvec mp = cross((rvec)(np2 - np1), np);
      mp = mp / sqrt(dot(mp, mp));

      // Line integral
      for (int j = 0; j < gausslegendre::N4; j++) {
        rvec sp = (np1 - 1e-6 * mp) * gausslegendre::x4[j][0] + 
                  (np2 - 1e-6 * mp) * gausslegendre::x4[j][1];
        double dLp = Lp * gausslegendre::w4[j];
        for (int q = 0; q < nQuadrature; ++q) {
          double dA = A * wQuadrature[q];
          // Convert source and observation coordinates to layered coordinates (for smooth evaluation)
          layeredCoords coordSmooth = 
                this->layeredUtils->cartesianToLayeredTab(svec[q], sp, interactionType);
          
          /* Evaluate smooth -- without quasistatic term -- tabulated Green's functions, 
          * namely: gTE, gTM, θz(gTE), θz(gTM), θzθz'(gTE), θzθz'(gTM), ...,
          * based on the interaction type (intra- or inter-layer).
          */
          auto gSmooth = (interactionType == "inter")
              ? this->evalGreenSmoothInter(coordSmooth.r, coordSmooth.Z1,
                                           coordSmooth.Z2, coeffIndex)
              : this->evalGreenSmoothIntra(interactionType, coordSmooth.r,
                                           coordSmooth.Z1, coordSmooth.Z2, coeffIndex);

          /* Evaluate quasistatic terms for Green's functions, namely: gTE_qs, 
          * gTM_qs, θz(gTE_qs), θzz(gTM_qs), ..., based on the interaction type
          * (intra- or inter-layer).
          */
          std::map<std::string, dcmplx> gQS;    // Quasistatic Green's functions
          std::map<std::string, dcmplx> gTotal; // Total Green's functions
          
          // Convert source and observation coordinates to layered coordinates (for quasistatic evaluation)
          layeredCoords coordQS = 
                this->layeredUtils->cartesianToLayeredQS(svec[q], sp, interactionType);
          // If no singularity cancellation is needed, use the tabulated quasistatic functions
          gQS = (interactionType == "inter")
              ? this->evalGreenQSInter(svec[q], sp, coordQS.r, coordQS.Z1, coordQS.R1)
              : this->evalGreenQSIntra(svec[q], interactionType, coordQS.r, coordQS.Z1, 
                                       coordQS.Z2, coordQS.R1, coordQS.R2);
          gTotal = gSmooth + gQS; // Green's function = smooth + quasistatic (+ operator 
                                  // is overloaded for this variable type)
          
          // Calculate the radial distance between source and observation points
          double rhoDist = sqrt((sp[0] - svec[q][0]) * (sp[0] - svec[q][0]) +
                                (sp[1] - svec[q][1]) * (sp[1] - svec[q][1]));
          
          if (rhoDist < 1e-10) rhoDist = 1e-10; // Avoid division by zero
          for (int k = 0; k < (int)fvec.size(); ++k) {
            RWGFun* f = fvec[k];
            rvec fs(f->Evaluate(svec[q]));
            double fKterm = (fs[0] * (sp[1] - svec[q][1]) - 
                             fs[1] * (sp[0] - svec[q][0])) / rhoDist;
            for (int l = 0; l < (int)fpvec.size(); ++l) {
              RWGFun* fp = fpvec[l];
              rvec fpsp(fp->Evaluate(sp));
              if (interactionType == "inter"){
                DK[k][l].second -=  // Opposite sign compared to Chew's formula 
                dLp * dA * dot(fpsp, mp) * fKterm * gTotal[TE + "rz2"];
              } else {
                DK[k][l].second -=  // Opposite sign compared to Chew's formula 
                dLp * dA * dot(fpsp, mp) * fKterm * gTotal[TE + "rz"];
              }
            }
          }
        }
      }
    }
  }
}

// Same as SameTriDK and additionally accumulates gradients.
void GreenFLayered3D::SameTriDKG(
    const std::vector<RWGFun*>& fvec, const std::vector<RWGFun*>& fpvec,
    std::vector<std::vector<std::pair<dcmplx, dcmplx>>>& DK,
    std::vector<std::vector<std::pair<cvec, cvec>>>& DKG) {
  // Since DKG should be called only for different domains, we assume no need of
  // singularity subtraction!

  Triangle *tPtr = fvec[0]->TrianglePtr(), *tpPtr = fpvec[0]->TrianglePtr();
  rvec r[3], rp[3];
  for (int i = 0; i < 3; i++) {
    r[i] = tPtr->Node(i);
    rp[i] = tpPtr->Node(i);
  }
  double A = tPtr->Area(), Ap = tpPtr->Area();

  int Ndunavant;
  double(*xdunavant)[3];
  double* wdunavant;
  if (NeedsAccurate) {
    Ndunavant = dunavant::N5;
    xdunavant = dunavant::x5;
    wdunavant = dunavant::w5;
  } else {
    Ndunavant = dunavant::N1;
    xdunavant = dunavant::x1;
    wdunavant = dunavant::w1;
  }

  std::vector<rvec> svec(Ndunavant), spvec(Ndunavant);
  for (int i = 0; i < Ndunavant; ++i) {
    svec[i] = r[0] * xdunavant[i][0] + r[1] * xdunavant[i][1] +
              r[2] * xdunavant[i][2];
    spvec[i] = rp[0] * xdunavant[i][0] + rp[1] * xdunavant[i][1] +
               rp[2] * xdunavant[i][2];
  }

  std::tuple<dcmplx, cvec, cdyad> allG;
  for (int i = 0; i < Ndunavant; ++i) {
    double dA = A * wdunavant[i];
    for (int j = 0; j < Ndunavant; ++j) {
      double dAp = Ap * wdunavant[j];
      allG = EvaluateAll(svec[i], spvec[j]);
      std::get<0>(allG) *= dA * dAp;
      std::get<1>(allG) *= dA * dAp;
      std::get<2>(allG) *= dA * dAp;
      for (int k = 0; k < (int)fvec.size(); ++k) {
        RWGFun* f = fvec[k];
        rvec fs(f->Evaluate(svec[i]));
        double divf(f->RWGPreFactor());
        for (int l = 0; l < (int)fpvec.size(); ++l) {
          RWGFun* fp = fpvec[l];
          rvec fpsp(fp->Evaluate(spvec[j]));
          double divfp(fp->RWGPreFactor());

          DK[k][l].first += std::get<0>(allG) *
                            (dot(fs, fpsp) - 1. / (k_B * k_B) * divf * divfp);
          DK[k][l].second += dot(fs, cross(std::get<1>(allG), (cvec)fpsp));
          DKG[k][l].first += std::get<1>(allG) *
                             (dot(fs, fpsp) - 1. / (k_B * k_B) * divf * divfp);
          DKG[k][l].second[0] +=
              dot(fs, cross(std::get<2>(allG)[0], (cvec)fpsp));
          DKG[k][l].second[1] +=
              dot(fs, cross(std::get<2>(allG)[1], (cvec)fpsp));
          DKG[k][l].second[2] +=
              dot(fs, cross(std::get<2>(allG)[2], (cvec)fpsp));
        }
      }
    }
  }
}

// Point–triangle delta/kappa contributions for RHS assembly.
void GreenFLayered3D::SameTriDeltaKappa(
    rvec r, const std::vector<RWGFun*>& fpvec,
    std::vector<std::pair<cvec, cvec>>& DeltaKappa)
{
  bool E_calc = true; /* Flag to check if DeltaKappaE or DeltaKappaH is calculated
                       * (according to DeltaKappa input initialization)
                       */
  auto &f = DeltaKappa[0].first;
  auto &s = DeltaKappa[0].second;
  bool firstIs100  = (f[0] == 1.0 && f[1] == 0.0 && f[2] == 0.0);
  bool firstIs000  = (f[0] == 0.0 && f[1] == 0.0 && f[2] == 0.0);
  bool secondIs000 = (s[0] == 0.0 && s[1] == 0.0 && s[2] == 0.0);
  // marker for "H‐only" initialization is first==cvec(1,0,0), second==cvec(0,0,0)
  if (firstIs100 && secondIs000) {
    // reset all entries to (0,0,0) so we can accumulate H‐terms
    for (auto &p : DeltaKappa) {
      p.first  = cvec(0.0, 0.0, 0.0);
      p.second = cvec(0.0, 0.0, 0.0);
    }
    E_calc = false;  // skip E, go straight to H
  }
  // if they initialized with anything else than (0,0,0)/(1,0,0), that's an error
  else if (!(firstIs000 && secondIs000)) 
  {
    throw std::runtime_error("DeltaKappa not initialized correctly.");
  }

  Triangle* Tp = fpvec[0]->TrianglePtr();
  rvec rp[3];
  for (int i = 0; i < 3; i++) {
    rp[i] = Tp->Node(i);
  }
  double Ap = Tp->Area();
  rvec Tc(Tp->Center());
  bool isOnTriangle = Tp->containsPoint(r);
  bool isOnInterface = false;
  for (int i = 0; i < zValsInterfaces.size(); i++) {
    if ( sqrt( ( r[2] - zValsInterfaces[i] ) * 
               ( r[2] - zValsInterfaces[i] ) ) < 1e-4 ) {
      isOnInterface = true; break;
    }
  }
  if (isOnTriangle && isOnInterface) {
    for (auto &p : DeltaKappa) {
      p.first  += cvec(0.0, 0.0, 0.0);
      p.second += cvec(0.0, 0.0, 0.0);
    }
    return;
  }
  bool needsSingCanc = !AboveThreshold(r, Tp);
  int nQuadrature;
  double(*xQuadrature)[3];
  double* wQuadrature;
  if (NeedsAccurate) {
    nQuadrature = needsSingCanc ? dunavant::N5 : dunavant::N1;
    xQuadrature = needsSingCanc ? dunavant::x5 : dunavant::x1;
    wQuadrature = needsSingCanc ? dunavant::w5 : dunavant::w1;
  } else {
    nQuadrature = dunavant::N1;
    xQuadrature = dunavant::x1;
    wQuadrature = dunavant::w1;
  }

  /* Set up the interaction type (intraTB, intraIn, inter),
   * coefficient and layer indices (nIndex, mIndex).
   * -----------------------------------------------------------
   * intraTB: intra-layer interaction in top/bottom layer,
   * intraIn: intra-layer interaction in the rest of the layers,
   * inter: inter-layer interactions,
   * nIndex: observation layer, mIndex: source layer
   */
  std::string interactionType = "unknown";
  int coeffIndex = -1;
  int nIndex = GetLayerIndex(r);
  int mIndex = GetLayerIndex(Tc);

  // Define the needed parameters for intra-layer and inter-layer interactions:
  // Source and observation layers are the same and are the top or bottom layer
  if ((nIndex == mIndex) && (nIndex == 0 || nIndex == epsilon.size() - 1)) {
    interactionType = "intraTB";
    auto it = std::find(intraTBIndex.begin(), intraTBIndex.end(), nIndex);
    if (it == intraTBIndex.end()) throw std::runtime_error("Invalid intraTB coeffIndex in SameTriDeltaKappa");
    coeffIndex = std::distance(intraTBIndex.begin(), it);
  }
  // Source and observation layers are the same but they are not the first or last layer
  else if ((nIndex == mIndex) && (nIndex != 0 && nIndex != epsilon.size() - 1)) {
    interactionType = "intraIn";
    auto it = std::find(intraPMIndex.begin(), intraPMIndex.end(), nIndex);
    if (it == intraPMIndex.end()) throw std::runtime_error("Invalid intraIn coeffIndex in SameTriDeltaKappa");
    coeffIndex = std::distance(intraPMIndex.begin(), it);
  }
  // Source and observation layers are different
  else if (nIndex != mIndex){
    interactionType = "inter";
    std::complex<int> target(nIndex, mIndex);
    auto it = std::find(interIndeces.begin(), interIndeces.end(), target);
    if (it == interIndeces.end()) throw std::runtime_error("Invalid inter coeffIndex in SameTriDeltaKappa");
    coeffIndex = std::distance(interIndeces.begin(), it);
  }
  else {
    std::cerr << "Something went wrong with the layers' indices during assembly evaluation" << std::endl;
  }

  dcmplx epsN(0.0, 0.0), epsM(0.0, 0.0), muN(0.0, 0.0), muM(0.0, 0.0);
  dcmplx kN(0.0, 0.0), kM(0.0, 0.0), kNM(0.0, 0.0), kMN(0.0, 0.0);
  std::string TE = "", TM = "";
  if (E_calc) { // Parameters for DK_E calculation
    epsN = epsilon[nIndex];  epsM = epsilon[mIndex];
    muN  = mu[nIndex];       muM  = mu[mIndex];
    kN   = k_L[nIndex];      kM   = k_L[mIndex];
    kNM  = 2.0 * PI * csqrt(epsN) * csqrt(muM) / wavelength;
    kMN  = 2.0 * PI * csqrt(epsM) * csqrt(muN) / wavelength;
    TE = "te"; TM = "tm";
  } else { // Parameters for DK_H calculation: duality principle (eps <-> mu, TE <-> TM)
    epsN = mu[nIndex];       epsM = mu[mIndex];
    muN  = epsilon[nIndex];  muM  = epsilon[mIndex];
    kN   = k_L[nIndex];      kM   = k_L[mIndex];
    kNM  = 2.0 * PI * csqrt(epsN) * csqrt(muM) / wavelength;
    kMN  = 2.0 * PI * csqrt(epsM) * csqrt(muN) / wavelength;
    TE = "tm"; TM = "te";
  }

  // Calculate reflected Green's function quasistatic part and its derivatives
  double dsp(1e-4);                                             // Displacement
  rvec rXplusDX(r[0] + dsp, r[1] + 0.0, r[2] + 0.0);            // Observation point + dx
  rvec rYplusDY(r[0] + 0.0, r[1] + dsp, r[2] + 0.0);            // Observation point + dy
  rvec rZplusDZ(r[0] + 0.0, r[1] + 0.0, r[2] + dsp);            // Observation point + dz
  std::vector<rvec> obsPoints = {rXplusDX, rYplusDY, rZplusDZ}; // Displaced observation points
  
  for (int i = 0; i < nQuadrature; ++i)
  {
    // quadrature point in the source triangle
    rvec spVec = rp[0] * xQuadrature[i][0] + rp[1] * xQuadrature[i][1] +
                 rp[2] * xQuadrature[i][2];

    double dAp = Ap * wQuadrature[i]; // weighted area of the observation triangle

    // Calculate the radial distance between source and observation points
    double rhoDist = sqrt((r[0] - spVec[0]) * (r[0] - spVec[0]) +
                          (r[1] - spVec[1]) * (r[1] - spVec[1]));
    if (rhoDist < 1e-10) rhoDist = 1e-10;

    // Calculate the unit vectors
    rvec unitRho = ( rvec( r[0] - spVec[0], r[1] - spVec[1], 0.0 ) ) / rhoDist;
    rvec unitZeta(0.0, 0.0, 1.0);
    rvec unitPhi = cross(unitZeta, unitRho);

    // Convert source and observation coordinates to layered coordinates (for smooth evaluation)
    layeredCoords coordSmooth = 
          this->layeredUtils->cartesianToLayeredTab(r, spVec, interactionType);
    
    /* Evaluate smooth -- without quasistatic term -- tabulated Green's functions, 
     * namely: gTE, gTM, θz(gTE), θz(gTM), θzθz'(gTE), θzθz'(gTM), ...,
     * based on the interaction type (intra- or inter-layer) and their spatial derivatives.
     */
    std::map<std::string, dcmplx> gSmooth;
    std::vector<std::map<std::string, dcmplx>> gSmoothGrad;
    gSmoothGrad.resize(obsPoints.size());
    gSmooth = (interactionType == "inter")
            ? this->evalGreenSmoothInter(coordSmooth.r, coordSmooth.Z1,
                                         coordSmooth.Z2, coeffIndex)
            : this->evalGreenSmoothIntra(interactionType, coordSmooth.r,
                                         coordSmooth.Z1, coordSmooth.Z2, coeffIndex);
    // derivatives calculation
    for (int i = 0; i < obsPoints.size(); ++i) {
      auto cSmoothtemp = this->layeredUtils->cartesianToLayeredTab(obsPoints[i], spVec,
                                                                   interactionType);
      auto gSmoothtemp = (interactionType == "inter")
                        ? this->evalGreenSmoothInter(cSmoothtemp.r, cSmoothtemp.Z1,
                                                     cSmoothtemp.Z2, coeffIndex)
                        : this->evalGreenSmoothIntra(interactionType, cSmoothtemp.r,
                                                     cSmoothtemp.Z1, cSmoothtemp.Z2, coeffIndex);
      gSmoothGrad[i] = (gSmoothtemp - gSmooth) / dsp;
    }

    /* Evaluate quasistatic terms for Green's functions, namely: gTE_qs, 
     * gTM_qs, θz(gTE_qs), θzz(gTM_qs), ..., based on the interaction type
     * (intra- or inter-layer) and their spatial derivatives.
     */
    std::map<std::string, dcmplx> gQS;                      // Quasistatic Green's functions
    std::vector<std::map<std::string, dcmplx>> gQSGrad;     // Quasistatic Green's functions gradient
    gQSGrad.resize(obsPoints.size());
    std::map<std::string, dcmplx> gTotal;                   // Total Green's functions
    std::vector<std::map<std::string, dcmplx>> gTotalGrad;  // Total Green's functions gradient
    gTotalGrad.resize(obsPoints.size());
    if (!needsSingCanc)
    {
      // Convert source and observation coordinates to layered coordinates (for quasistatic evaluation)
      auto coordQS = 
           this->layeredUtils->cartesianToLayeredQS(r, spVec, interactionType);
      
      // If no singularity cancellation is needed, use the typical quasistatic functions
      gQS = (interactionType == "inter")
          ? this->evalGreenQSInter(r, spVec, coordQS.r, coordQS.Z1, coordQS.R1)
          : this->evalGreenQSIntra(r, interactionType, coordQS.r, coordQS.Z1, 
                                   coordQS.Z2, coordQS.R1, coordQS.R2);
      
      gTotal = gSmooth + gQS; /* Green's function = smooth + quasistatic (+ operator 
                               * is overloaded for this variable type) */
      
      // derivatives calculation
      for (int i = 0; i < obsPoints.size(); ++i) {
        auto cQStemp = this->layeredUtils->cartesianToLayeredQS(obsPoints[i], spVec, interactionType);
        auto gQStemp = (interactionType=="inter")
                     ? evalGreenQSInter(obsPoints[i], spVec, cQStemp.r, cQStemp.Z1, cQStemp.R1)
                     : evalGreenQSIntra(obsPoints[i], interactionType, cQStemp.r, cQStemp.Z1,
                                        cQStemp.Z2, cQStemp.R1, cQStemp.R2);
        gQSGrad[i] = (gQStemp - gQS) / dsp;
        gTotalGrad[i] = gSmoothGrad[i] + gQSGrad[i];
      }
    } else {
      // If singularity cancellation is needed, use the analytical quasistatic functions
      gTotal = gSmooth; // Only smooth Green's functions are used
      for (int i = 0; i < obsPoints.size(); ++i) {
        gTotalGrad[i] = gSmoothGrad[i];
      }
    }

    // build a cvec of the three partials
    auto gradGTotal = [&]( const std::string &key ) -> cvec {
      return cvec{
        gTotalGrad[0].at(key),
        gTotalGrad[1].at(key),
        gTotalGrad[2].at(key)
      };
    };

    for (int j = 0; j < (int)fpvec.size(); ++j) {
      rvec fpTot    = fpvec[j]->Evaluate(spVec);
      double fpRho  = dot(fpTot, unitRho);
      double fpPhi  = dot(fpTot, unitPhi);
      double fpZeta = dot(fpTot, unitZeta);
      double divfp  = fpvec[j]->RWGPreFactor();

      // compute DeltaKappa
      DeltaKappa[j].first +=
      dAp * ( (-gradGTotal(TM + "zz") / (kNM * kNM) + gradGTotal(TE) ) * divfp  +
              ( gTotal[TM] * kMN * kMN - gTotal[TE + "zz"]) * unitZeta * fpZeta +
              (fpRho * unitRho + fpPhi * unitPhi) * gTotal[TE + "s"] );
      DeltaKappa[j].second +=
      dAp * ( gTotal[TE + "r"] * kM * kM * unitPhi *  fpZeta -
              gTotal[TM + "r"] * kMN * kMN * unitZeta * fpPhi );
      if (interactionType == "inter") {
        DeltaKappa[j].first -= 
        dAp * ( ( gTotal[TM + "z2"] * (muN / muM) + gTotal[TE + "z1"] ) * unitZeta * divfp -
                ( gradGTotal(TM + "z1") * (epsM / epsN) + gradGTotal(TE + "z2") ) * fpZeta );
        DeltaKappa[j].second +=
        dAp * ( gTotal[TM + "rz1"] * (epsM / epsN) * unitPhi * fpRho / rhoDist -
                gradGTotal(TM + "rz1") * (epsM / epsN) * fpPhi -
                gTotal[TE + "rz2"] * unitPhi * divfp );
      } else {
        DeltaKappa[j].first -= 
        dAp * ( ( gTotal[TM + "z"] * (muN / muM) + gTotal[TE + "z"] ) * unitZeta * divfp -
                ( gradGTotal(TM + "z") * (epsM / epsN) + gradGTotal(TE + "z") ) * fpZeta );
        DeltaKappa[j].second +=
        dAp * ( gTotal[TM + "rz"] * (epsM / epsN) * unitPhi * fpRho / rhoDist -
                gradGTotal(TM + "rz") * (epsM / epsN) * fpPhi -
                gTotal[TE + "rz"] * unitPhi * divfp );
      }
    }
  }
  if (needsSingCanc) {
    TriQuadPol triQuad(
      Tp,                 // pointer to source triangle
      r,                  // observation point
      gausslegendre::N3,  // number of Gauss-Legendre points for the 1D rule
      gausslegendre::x3,  // table of 1D Gauss-Legendre points
      gausslegendre::w3   // table of 1D Gauss-Legendre weights
    );
    // Get the quadrature points and weights
    std::vector<rvec> quadPoints = triQuad.getQuadPoints();
    std::vector<double> quadWeights = triQuad.getQuadWeights();
#ifdef DEBUG
    writeVectorToFile(quadPoints, quadWeights,
                     "quadPoints_rho_" + std::to_string(r[0]) + "_" + 
                                         std::to_string(r[1]) + "_" + 
                                         std::to_string(r[2]) + 
                      "_Tp_" + std::to_string((uintptr_t)Tp) + ".txt");
#endif
    // Loop over the quadrature points
    for (int idx = 0; idx < quadWeights.size(); ++idx)
    {
      // quadrature point in the source triangle
      rvec spVec = rp[0] * quadPoints[idx][0] + rp[1] * quadPoints[idx][1] + 
                   rp[2] * quadPoints[idx][2];
      double dW = quadWeights[idx];  // weights for the observation triangle
      if (dW < 1e-10) continue;
      // Calculate the radial distance between source and observation points
      double rhoDist = sqrt((r[0] - spVec[0]) * (r[0] - spVec[0]) +
                            (r[1] - spVec[1]) * (r[1] - spVec[1]));
      if (rhoDist < 1e-10) rhoDist = 1e-10;

      // Calculate the unit vectors
      rvec unitRho = ( rvec( r[0] - spVec[0], r[1] - spVec[1], 0.0 ) ) / rhoDist;
      rvec unitZeta(0.0, 0.0, 1.0);
      rvec unitPhi = cross(unitZeta, unitRho);

      /* Evaluate quasistatic terms for Green's functions, namely: gTE_qs, 
       * gTM_qs, θz(gTE_qs), θzz(gTM_qs), ..., based on the interaction type
       * (intra- or inter-layer) and their spatial derivatives.
       */
      std::map<std::string, dcmplx> gQS;                  // Quasistatic Green's functions
      std::vector<std::map<std::string, dcmplx>> gQSGrad; // Quasistatic Green's functions gradient
      gQSGrad.resize(obsPoints.size());

      // Convert source and observation coordinates to layered coordinates (for quasistatic evaluation)
      auto coordQS = layeredUtils->cartesianToLayeredQS(r, spVec, interactionType);

      // Since singularity cancellation is needed, use the modified quasistatic functions
      gQS = (interactionType=="inter")
          ? evalGreenQSInter(r, spVec, coordQS.r, coordQS.Z1, coordQS.R1, needsSingCanc)
          : evalGreenQSIntra(r, interactionType, coordQS.r, coordQS.Z1,
                             coordQS.Z2, coordQS.R1, coordQS.R2, needsSingCanc);
      
      // derivatives calculation
      for (int i = 0; i < obsPoints.size(); ++i) {
        auto cQStemp = layeredUtils->cartesianToLayeredQS(obsPoints[i], spVec, interactionType);
        auto gQStemp = (interactionType=="inter")
                     ? evalGreenQSInter(obsPoints[i], spVec, cQStemp.r, cQStemp.Z1,
                                        cQStemp.R1, needsSingCanc)
                     : evalGreenQSIntra(obsPoints[i], interactionType, cQStemp.r, cQStemp.Z1,
                                        cQStemp.Z2, cQStemp.R1, cQStemp.R2, needsSingCanc);
        gQSGrad[i] = (gQStemp - gQS) / dsp;
      }

      // build a cvec of the three partials
      auto gradGQS = [&]( const std::string &key ) -> cvec {
        return cvec{
          gQSGrad[0].at(key),
          gQSGrad[1].at(key),
          gQSGrad[2].at(key)
        };
      };

      for (int j = 0; j < (int)fpvec.size(); ++j) {
        rvec fpTot    = fpvec[j]->Evaluate(spVec);
        double fpRho  = dot(fpTot, unitRho);
        double fpPhi  = dot(fpTot, unitPhi);
        double fpZeta = dot(fpTot, unitZeta);
        double divfp  = fpvec[j]->RWGPreFactor();
        
        // compute DeltaKappa
        if (interactionType == "intraTB") {
          DeltaKappa[j].first += 
          dW * ( -( gQS[TM + "Refl"] + gQS[TE + "Refl"] ) * 
                  ( gQS["z"] * unitZeta * divfp - gradGQS("z") * fpZeta ) -
                  gQS[TM + "Refl"] * gradGQS("zz") * divfp / (kNM * kNM) -
                  gQS[TE + "Refl"] * gQS["zz"] * unitZeta * fpZeta + 
                  gQS[TE + "Refl"] * gQS["s"] * ( unitRho * fpRho + unitPhi * fpPhi ) );
          DeltaKappa[j].second += 
          dW * ( gQS[TM + "Refl"] * 
                 ( gQS["rz"] * unitPhi * fpRho / rhoDist - gradGQS("rz") * fpPhi ) -
                 gQS[TM + "Refl"] * gQS["r"] * unitZeta * fpPhi * kMN * kMN -
                 gQS[TE + "Refl"] * gQS["rz"] * unitPhi * divfp + 
                 gQS[TE + "Refl"] * gQS["r"] * unitPhi * fpZeta * kM * kM );
        }
        else if (interactionType == "intraIn") {
          for (int k = 0; k < 2; ++k) {
            std::string intNum = std::to_string(k + 1); // lower (= 1) or upper (= 2) layer contribution
            DeltaKappa[j].first += 
            dW * ( -( gQS[TM + "Refl" + intNum] + gQS[TE + "Refl" + intNum] ) * 
                    ( gQS["z" + intNum] * unitZeta * divfp - gradGQS("z" + intNum) * fpZeta ) -
                    gQS[TM + "Refl" + intNum] * gradGQS("zz" + intNum) * divfp / (kNM * kNM) -
                    gQS[TE + "Refl" + intNum] * gQS["zz" + intNum] * unitZeta * fpZeta + 
                    gQS[TE + "Refl" + intNum] * 
                    gQS["s" + intNum] * ( unitRho * fpRho + unitPhi * fpPhi ) );
            DeltaKappa[j].second += 
            dW * ( gQS[TM + "Refl" + intNum] * 
                   ( gQS["rz" + intNum] * unitPhi * fpRho / rhoDist - gradGQS("rz" + intNum) * fpPhi ) -
                   gQS[TM + "Refl" + intNum] * gQS["r" + intNum] * unitZeta * fpPhi * kMN * kMN -
                   gQS[TE + "Refl" + intNum] * gQS["rz" + intNum] * unitPhi * divfp + 
                   gQS[TE + "Refl" + intNum] * gQS["r" + intNum] * unitPhi * fpZeta * kM * kM );
          }
        }
        else if (interactionType == "inter") {
          DeltaKappa[j].first +=
          dW * ( - gQS[TM + "Trans"] * gradGQS("zz") * divfp / (kNM * kNM) - 
                   gQS[TM + "Trans"] * gQS["z2"] * unitZeta * divfp * (muN / muM) - 
                   gQS[TE + "Trans"] * gQS["z1"] * unitZeta * divfp + 
                   gQS[TM + "Trans"] * gradGQS("z1") * fpZeta * (epsM / epsN) + 
                   gQS[TE + "Trans"] * gradGQS("z2") * fpZeta - 
                   gQS[TE + "Trans"] * gQS["zz"] * unitZeta * fpZeta + 
                   gQS[TE + "Trans"] * gQS["s"] * ( unitRho * fpRho + unitPhi * fpPhi ) );
          DeltaKappa[j].second +=
          dW * ( gQS[TM + "Trans"] * 
                 ( ( gQS["rz1"] * unitPhi * fpRho / rhoDist - gradGQS("rz1") * fpPhi ) * (epsM / epsN) - 
                   gQS["r"] * unitZeta * fpPhi * kMN * kMN ) - 
                 gQS[TE + "Trans"] * 
                 ( gQS["rz2"] * unitPhi * divfp - gQS["r"] * unitPhi * fpZeta * kM * kM) );
        }
        else {
          throw std::runtime_error("Unknown interaction type: " + interactionType);
        }
      }
    }
  }
}

// Computes Fresnel reflection and transmission coefficients.
std::vector<dcmplx> GreenFLayered3D::FresnelCoeff(int layerIndexI, int layerIndexJ, 
                                                  dcmplx krho, const std::string& waves) const {
  return this->layeredUtils->FresnelCoeff(layerIndexI, layerIndexJ, krho, waves);
}

// Fills the Green's function table for intra-layer interactions (top/bottom layers).
void GreenFLayered3D::FillIntraTopBottom(const dcmplx& rq_te, const dcmplx& rq_tm, const Grid& inGrid) {
  std::vector<std::vector<double>> rtab, ztab;      // 2D grids
  size_t nr = inGrid.r.size();                      // Number of r values
  size_t nz = inGrid.z1.size();                     // Number of z values

  rtab.resize(nr, std::vector<double>(nz)); // Initialize r grid
  ztab.resize(nr, std::vector<double>(nz)); // Initialize z grid

  // Generate grids
  for (size_t i = 0; i < nr; ++i) {
    for (size_t j = 0; j < nz; ++j) {
      rtab[i][j] = inGrid.r.at(i);    // Each column has rtab values
      ztab[i][j] = inGrid.z1.at(j);   // Each row has ztab values
    }
  }

  // Evaluate Sommerfeld integral (returns interpolation coefficients & properties)
  sommIntIntraTB.push_back(
      new SommerfeldIntegrator(k_L, inGrid.i1, epsilon, mu, zValsInterfaces, thickness));
  coeffIntraTB.push_back( 
      sommIntIntraTB.back()->Evaluate({rq_te}, {rq_tm}, rtab, ztab, true));

  // Save material index
  intraTBIndex.push_back(inGrid.i1);
}

// Fills the Green's function table for intra-layer interactions (internal layers).
void GreenFLayered3D::FillIntraInside(const std::vector<dcmplx>& rq_te, 
                                      const std::vector<dcmplx>& rq_tm, 
                                      const Grid& inGrid) {
  std::vector<std::vector<double>> rtab, z1tab, z2tab;  // 2D grids
  size_t nr = inGrid.r.size();                          // Number of r values
  size_t nz = inGrid.z1.size();                         // Number of z1-z2 values

  rtab.resize(nz, std::vector<double>(nr));  // Initialize r grid
  z1tab.resize(nz, std::vector<double>(nr)); // Initialize z1 grid
  z2tab.resize(nz, std::vector<double>(nr)); // Initialize z2 grid

  // Generate grids
  for (size_t i = 0; i < nr; ++i) {
    for (size_t j = 0; j < nz; ++j) {
      rtab[i][j] = inGrid.r.at(i);    // Each column has rtab values
      z1tab[i][j] = inGrid.z1.at(j);  // Each row has z1tab values
      z2tab[i][j] = inGrid.z2.at(j);  // Each row has z2tab values
    }
  }
  // Evaluate Sommerfeld integral (returns interpolation coefficients & properties)
  sommIntIntraMinus.push_back(
      new SommerfeldIntegrator(k_L, inGrid.i1, epsilon, mu, zValsInterfaces, thickness));
  sommIntIntraPlus.push_back(
      new SommerfeldIntegrator(k_L, inGrid.i1, epsilon, mu, zValsInterfaces, thickness));
  coeffIntraMinus.push_back(
      sommIntIntraMinus.back()->Evaluate(rq_te, rq_tm, rtab, z1tab, true, "minus"));
  coeffIntraPlus.push_back(
      sommIntIntraPlus.back()->Evaluate(rq_te, rq_tm, rtab, z2tab, true, "plus"));

  // Save material index
  intraPMIndex.push_back(inGrid.i1);
}

// Fills the Green's function table for inter-layer interactions.
void GreenFLayered3D::FillInter(const dcmplx& tq_te, const dcmplx& tq_tm, 
                                const Grid& inGrid) {
  std::vector<std::vector<std::vector<double>>> rtab, ztab1, ztab2; // 3D grids
  size_t nr  = inGrid.r.size();                                     // Number of r values
  size_t nz1 = inGrid.z1.size();                                    // Number of z1 values
  size_t nz2 = inGrid.z2.size();                                    // Number of z2 values
  
  // Resize to 3D shape (nz2, nz1, nr)
  rtab.resize(nz2, std::vector<std::vector<double>>(nz1, std::vector<double>(nr)));
  ztab1.resize(nz2, std::vector<std::vector<double>>(nz1, std::vector<double>(nr)));
  ztab2.resize(nz2, std::vector<std::vector<double>>(nz1, std::vector<double>(nr)));

  // Generate 3D grids
  for (size_t i = 0; i < nr; ++i) {       // Iterate over r dimension
    for (size_t j = 0; j < nz1; ++j) {    // Iterate over z1 dimension
      for (size_t k = 0; k < nz2; ++k) {  // Iterate over z2 dimension
        rtab[i][j][k] = inGrid.r.at(i);   // rtab values in the last dimension
        ztab1[i][j][k] = inGrid.z1.at(j); // ztab1 values in the second dimension
        ztab2[i][j][k] = inGrid.z2.at(k); // ztab2 values in the first dimension
      }
    }
  }

  // Evaluate Sommerfeld integral (returns interpolation coefficients & properties)
  sommIntInter.push_back(
      new SommerfeldIntegrator(k_L, inGrid.i1 , inGrid.i2, epsilon, 
                               mu, zValsInterfaces, thickness));
  coeffInter.push_back(
      sommIntInter.back()->Evaluate(tq_te, tq_tm, rtab, ztab1, ztab2, true));

  // Save material index
  interIndeces.push_back( std::complex<int>(inGrid.i1, inGrid.i2) );
}


// Fills all interpolation tables with Green's function values.
void GreenFLayered3D::FillGreenFTable() {
  // Iterate over all tabulation grids
  for (int i = 0; i < tabulationGrids.size(); i++) 
  {
    std::vector<dcmplx> fresnelQS;  // Fresnel reflection and transmission coefficients (quasistatic)

    // Source and observation layers are the same and are the top or bottom layer
    if (tabulationGrids[i].i1 == tabulationGrids[i].i2 &&
       (tabulationGrids[i].i1 == 0 || tabulationGrids[i].i1 == epsilon.size() - 1)) {        
      if (tabulationGrids[i].i1 == 0) { // Bottom layer
        std::vector<dcmplx> temp = FresnelCoeff(0, 1, dcmplx(10000000000.0, 0.0), "TE");
        fresnelQS.push_back(temp[0]);   // Reflection coefficient for TE wave
        temp = FresnelCoeff(0, 1, dcmplx(10000000000.0, 0.0), "TM");
        fresnelQS.push_back(temp[0]);   // Reflection coefficient for TM wave
      }
      else { // Top layer
        std::vector<dcmplx> temp = FresnelCoeff(epsilon.size() - 1, epsilon.size() - 2, 
                                                dcmplx(10000000000.0, 0.0), "TE");
        fresnelQS.push_back(temp[0]); // Reflection coefficient for TE wave
        temp = FresnelCoeff(epsilon.size() - 1, epsilon.size() - 2, 
                            dcmplx(10000000000.0, 0.0), "TM");
        fresnelQS.push_back(temp[0]); // Reflection coefficient for TM wave
      }

      // Ensure that monotonic axes are ascending
      bool revR = tabulationGrids[i].r.front()  > tabulationGrids[i].r.back();
      bool revZ = tabulationGrids[i].z1.front() > tabulationGrids[i].z1.back();
      if (revR) std::reverse(tabulationGrids[i].r.begin(),  tabulationGrids[i].r.end());
      if (revZ) std::reverse(tabulationGrids[i].z1.begin(), tabulationGrids[i].z1.end());

      // Fill the Green's function table for the top and bottom layers
      FillIntraTopBottom(fresnelQS[0], fresnelQS[1], tabulationGrids[i]);
    }
    // Source and observation layers are the same but they are not the first or last layer
    else if (tabulationGrids[i].i1 == tabulationGrids[i].i2) {
      std::vector<dcmplx> temp = FresnelCoeff(tabulationGrids[i].i1, tabulationGrids[i].i1 - 1, 
                                              dcmplx(10000000000.0, 0.0), "TE");
      fresnelQS.push_back(temp[0]); // Reflection coefficient R_{i1, i1 - 1} for TE wave
      temp = FresnelCoeff(tabulationGrids[i].i1, tabulationGrids[i].i1 + 1, 
                          dcmplx(10000000000.0, 0.0), "TE");
      fresnelQS.push_back(temp[0]); // Reflection coefficient R_{i1, i1 + 1} for TE wave
      temp = FresnelCoeff(tabulationGrids[i].i1, tabulationGrids[i].i1 - 1, 
                          dcmplx(10000000000.0, 0.0), "TM");
      fresnelQS.push_back(temp[0]); // Reflection coefficient R_{i1, i1 - 1} for TM wave
      temp = FresnelCoeff(tabulationGrids[i].i1, tabulationGrids[i].i1 + 1, 
                          dcmplx(10000000000.0, 0.0), "TM");
      fresnelQS.push_back(temp[0]); // Reflection coefficient R_{i1, i1 + 1} for TM wave

      // Ensure that monotonic axes are ascending
      bool revR  = tabulationGrids[i].r.front()  > tabulationGrids[i].r.back();
      bool revZ1 = tabulationGrids[i].z1.front() > tabulationGrids[i].z1.back();
      bool revZ2 = tabulationGrids[i].z2.front() > tabulationGrids[i].z2.back();
      if (revR)  std::reverse(tabulationGrids[i].r.begin(),  tabulationGrids[i].r.end());
      if (revZ1) std::reverse(tabulationGrids[i].z1.begin(), tabulationGrids[i].z1.end());
      if (revZ2) std::reverse(tabulationGrids[i].z2.begin(), tabulationGrids[i].z2.end());

      // Fill the Green's function table for the inside layers
      FillIntraInside({fresnelQS[0], fresnelQS[1]}, {fresnelQS[2], fresnelQS[3]}, tabulationGrids[i]);
    }
    // Source and observation layers are different
    else if (tabulationGrids[i].i1 != tabulationGrids[i].i2){
      if (std::abs(tabulationGrids[i].i1 - tabulationGrids[i].i2) == 1) {
        std::vector<dcmplx> temp = FresnelCoeff(tabulationGrids[i].i2, tabulationGrids[i].i1, 
                                                dcmplx(10000000000.0, 0.0), "TE");
        fresnelQS.push_back(temp[1]); // Transmission coefficient for TE wave
        temp = FresnelCoeff(tabulationGrids[i].i2, tabulationGrids[i].i1, 
                            dcmplx(10000000000.0, 0.0), "TM");
        fresnelQS.push_back(temp[1]); // Transmission coefficient for TM wave
      }
      else {
        fresnelQS.push_back(dcmplx(0.0, 0.0)); // Transmission coefficient for TE wave
        fresnelQS.push_back(dcmplx(0.0, 0.0)); // Transmission coefficient for TM wave
      }

      // Ensure that monotonic axes are ascending
      bool revR  = tabulationGrids[i].r.front()  > tabulationGrids[i].r.back();
      bool revZ1 = tabulationGrids[i].z1.front() > tabulationGrids[i].z1.back();
      bool revZ2 = tabulationGrids[i].z2.front() > tabulationGrids[i].z2.back();
      if (revR)  std::reverse(tabulationGrids[i].r.begin(),  tabulationGrids[i].r.end());
      if (revZ1) std::reverse(tabulationGrids[i].z1.begin(), tabulationGrids[i].z1.end());
      if (revZ2) std::reverse(tabulationGrids[i].z2.begin(), tabulationGrids[i].z2.end());

      // Fill the Green's function table for the inter-layer interactions
      FillInter(fresnelQS[0], fresnelQS[1], tabulationGrids[i]);
    }
    else {
      std::cerr << "Somenthing went wrong with the layers' indices during tabulation" << std::endl;
    }
  } 
}

// Evaluate smooth Green's function components for **intra-layer** interactions.
std::map<std::string, dcmplx> GreenFLayered3D::evalGreenSmoothIntra(const std::string& interactionType, double R, 
                                                                    double Z1, double Z2, int coeffIndex) {
    std::map<std::string,dcmplx> out; // Prepare the output map

    // all of the intra‐map keys:
    static const std::vector<std::string> keys = {"te","tm","tez","tmz","tezz","tmzz",
                                                  "tes","tms","ter","tmr","terz","tmrz"};

    if (interactionType == "intraTB") {
        // top/bottom: single table at (R,Z1)
        for (auto const& key : keys) {
            auto& tpl = coeffSelectorIntra.at(key)(coeffIntraTB.at(coeffIndex));
            double re = evaluateSinglePoint(R, Z1, tpl, "real");
            double im = evaluateSinglePoint(R, Z1, tpl, "imag");
            out[key] = dcmplx(re, im);
        }
    }
    else { // intraIn: minus at Z1, plus at Z2
        for (auto const& key : keys) {
            // minus
            auto& tplMinus = coeffSelectorIntra.at(key)(coeffIntraMinus.at(coeffIndex));
            double reMinus = evaluateSinglePoint(R, Z1, tplMinus, "real");
            double imMinus = evaluateSinglePoint(R, Z1, tplMinus, "imag");
            
            // plus
            auto& tplPlus = coeffSelectorIntra.at(key)(coeffIntraPlus.at(coeffIndex));
            double rePlus = evaluateSinglePoint(R, Z2, tplPlus, "real");
            double imPlus = evaluateSinglePoint(R, Z2, tplPlus, "imag");
            out[key] = dcmplx(reMinus + rePlus, imMinus + imPlus); // Combine minus and plus
        }
    }

    return out;
}

// Evaluate smooth Green's function components for **inter-layer** interactions.
std::map<std::string, dcmplx> GreenFLayered3D::evalGreenSmoothInter(double R, double Z1, double Z2, 
                                                                   int coeffIndex) {    
    std::map<std::string,dcmplx> out; // Prepare the output map

    // Capture the raw, unscaled values once:
    double rho = R  * coeffInter.at(coeffIndex).scales[0];
    double z1  = Z1 * coeffInter.at(coeffIndex).scales[1];
    double z2  = Z2 * coeffInter.at(coeffIndex).scales[2];
    
    // all of the inter‐map keys:
    static const std::vector<std::string> keys = {"te","tm","tez1","tez2","tmz1","tmz2",
                                                  "tezz","tmzz","tes","tms","ter","tmr",
                                                  "terz1","terz2","tmrz1","tmrz2"};
    for (auto const& key : keys) {
        auto& tpl = coeffSelectorInter.at(key)(coeffInter.at(coeffIndex));
        double re = this->evaluateSinglePoint(rho, z1, z2, tpl, "real");
        double im = this->evaluateSinglePoint(rho, z1, z2, tpl, "imag");
        out[key] = dcmplx(re, im);
    }

    return out;
}

// Evaluate quasistatic Green's function components for **intra-layer** interactions.
std::map<std::string, dcmplx> GreenFLayered3D::evalGreenQSIntra(
  const rvec pos, const std::string& interactionType, double rho,
  double Z1, double Z2, double R1, double R2, bool singCancFields)
{
  std::map<std::string,dcmplx> out; // Prepare the output map
  std::vector<dcmplx> fresnelQS;    // Quasistatic Fresnel reflection coefficients

  if (interactionType == "intraTB") { // top/bottom: single table at (R,Z1)
    double dir = 0.0; // Propagation direction of reflected waves
    if (this->GetLayerIndex(pos) == 0) {
        dir = -1.0; // Bottom layer
        std::vector<dcmplx> temp = FresnelCoeff(0, 1, dcmplx(10000000000.0, 0.0), "TE");
        fresnelQS.push_back(temp[0]);   // Reflection coefficient for TE wave
        temp = FresnelCoeff(0, 1, dcmplx(10000000000.0, 0.0), "TM");
        fresnelQS.push_back(temp[0]);   // Reflection coefficient for TM wave
    } else {
        dir = 1.0;  // Top layer
        std::vector<dcmplx> temp = FresnelCoeff(epsilon.size() - 1, epsilon.size() - 2, 
                                                dcmplx(10000000000.0, 0.0), "TE");
        fresnelQS.push_back(temp[0]); // Reflection coefficient for TE wave
        temp = FresnelCoeff(epsilon.size() - 1, epsilon.size() - 2, 
                            dcmplx(10000000000.0, 0.0), "TM");
        fresnelQS.push_back(temp[0]); // Reflection coefficient for TM wave
    }

    // quasistatic Green's functions terms
    dcmplx termZ = dir * log(Z1 + R1) / (4.0 * PI);
    dcmplx termZZ_S = 1 / (4.0 * PI * R1);
    dcmplx termR = - (R1 - Z1) / (4.0 * PI * rho);
    dcmplx termRZ = dir * (1.0 - Z1 / R1) / (4.0 * PI * rho);

    // form depending on singularity cancellation for E, H fields calculation
    if (singCancFields) {
      out["z"]  = termZ;
      out["zz"] = termZZ_S;
      out["s"]  = termZZ_S;
      out["r"]  = termR;
      out["rz"] = termRZ;
      out["teRefl"] = fresnelQS[0];
      out["tmRefl"] = fresnelQS[1];
    } else {
      out["tez"]  = (std::abs(fresnelQS[0]) > 1.0e-10) ? fresnelQS[0] * termZ : 0.0;
      out["tmz"]  = (std::abs(fresnelQS[1]) > 1.0e-10) ? fresnelQS[1] * termZ : 0.0;
      out["tezz"] = (std::abs(fresnelQS[0]) > 1.0e-10) ? fresnelQS[0] * termZZ_S : 0.0;
      out["tmzz"] = (std::abs(fresnelQS[1]) > 1.0e-10) ? fresnelQS[1] * termZZ_S : 0.0;
      out["tes"]  = (std::abs(fresnelQS[0]) > 1.0e-10) ? fresnelQS[0] * termZZ_S : 0.0;
      out["tms"]  = (std::abs(fresnelQS[1]) > 1.0e-10) ? fresnelQS[1] * termZZ_S : 0.0;
      out["ter"]  = (std::abs(fresnelQS[0]) > 1.0e-10) ? fresnelQS[0] * termR : 0.0;
      out["tmr"]  = (std::abs(fresnelQS[1]) > 1.0e-10) ? fresnelQS[1] * termR : 0.0;
      out["terz"] = (std::abs(fresnelQS[0]) > 1.0e-10) ? fresnelQS[0] * termRZ : 0.0;
      out["tmrz"] = (std::abs(fresnelQS[1]) > 1.0e-10) ? fresnelQS[1] * termRZ : 0.0;
    }
  }
  else { // intraIn: minus at Z1, plus at Z2
    std::vector<dcmplx> temp = FresnelCoeff(this->GetLayerIndex(pos), this->GetLayerIndex(pos) - 1, 
                                            dcmplx(10000000000.0, 0.0), "TE");
    fresnelQS.push_back(temp[0]); // Reflection coefficient R_{i1, i1 - 1} for TE wave
    temp = FresnelCoeff(this->GetLayerIndex(pos), this->GetLayerIndex(pos) + 1, 
                        dcmplx(10000000000.0, 0.0), "TE");
    fresnelQS.push_back(temp[0]); // Reflection coefficient R_{i1, i1 + 1} for TE wave
    temp = FresnelCoeff(this->GetLayerIndex(pos), this->GetLayerIndex(pos) - 1, 
                        dcmplx(10000000000.0, 0.0), "TM");
    fresnelQS.push_back(temp[0]); // Reflection coefficient R_{i1, i1 - 1} for TM wave
    temp = FresnelCoeff(this->GetLayerIndex(pos), this->GetLayerIndex(pos) + 1, 
                        dcmplx(10000000000.0, 0.0), "TM");
    fresnelQS.push_back(temp[0]); // Reflection coefficient R_{i1, i1 + 1} for TM wave
    
    // quasistatic Green's functions terms
    dcmplx termZ1    =  1.0 * log(Z1 + R1) / (4.0 * PI);
    dcmplx termZ2    = -1.0 * log(Z2 + R2) / (4.0 * PI);
    dcmplx termZZ_S1 =  1.0 / (4.0 * PI * R1);
    dcmplx termZZ_S2 =  1.0 / (4.0 * PI * R2);
    dcmplx termR1    = -(R1 - Z1) / (4.0 * PI * rho);
    dcmplx termR2    = -(R2 - Z2) / (4.0 * PI * rho);
    dcmplx termRZ1   =  (1.0 - Z1 / R1) / (4.0 * PI * rho);
    dcmplx termRZ2   = -(1.0 - Z2 / R2) / (4.0 * PI * rho);

    // form depending on singularity cancellation for E, H fields calculation
    if (singCancFields) {
      out["z1"]  = termZ1;              out["z2"]  = termZ2;
      out["zz1"] = termZZ_S1;           out["zz2"] = termZZ_S2;
      out["s1"]  = termZZ_S1;           out["s2"]  = termZZ_S2;
      out["r1"]  = termR1;              out["r2"]  = termR2;
      out["rz1"] = termRZ1;             out["rz2"] = termRZ2;
      out["teRefl1"] = fresnelQS[0];    out["teRefl2"] = fresnelQS[1];
      out["tmRefl1"] = fresnelQS[2];    out["tmRefl2"] = fresnelQS[3];
    } else {
      out["tez"] = (std::abs(fresnelQS[0]) + std::abs(fresnelQS[1]) > 1.0e-10) ? 
                    fresnelQS[0] * termZ1 + fresnelQS[1] * termZ2 : 0.0;
      out["tmz"] = (std::abs(fresnelQS[2]) + std::abs(fresnelQS[3]) > 1.0e-10) ? 
                    fresnelQS[2] * termZ1 + fresnelQS[3] * termZ2 : 0.0;
      out["tezz"] = (std::abs(fresnelQS[0]) + std::abs(fresnelQS[1]) > 1.0e-10) ? 
                    fresnelQS[0] * termZZ_S1 + fresnelQS[1] * termZZ_S2 : 0.0;
      out["tmzz"] = (std::abs(fresnelQS[2]) + std::abs(fresnelQS[3]) > 1.0e-10) ? 
                    fresnelQS[2] * termZZ_S1 + fresnelQS[3] * termZZ_S2 : 0.0;
      out["tes"]  = (std::abs(fresnelQS[0]) + std::abs(fresnelQS[1]) > 1.0e-10) ? 
                    fresnelQS[0] * termZZ_S1 + fresnelQS[1] * termZZ_S2 : 0.0;
      out["tms"]  = (std::abs(fresnelQS[2]) + std::abs(fresnelQS[3]) > 1.0e-10) ? 
                    fresnelQS[2] * termZZ_S1 + fresnelQS[3] * termZZ_S2 : 0.0;
      out["ter"] = (std::abs(fresnelQS[0]) + std::abs(fresnelQS[1]) > 1.0e-10) ? 
                    fresnelQS[0] * termR1 + fresnelQS[1] * termR2 : 0.0;
      out["tmr"] = (std::abs(fresnelQS[2]) + std::abs(fresnelQS[3]) > 1.0e-10) ? 
                    fresnelQS[2] * termR1 + fresnelQS[3] * termR2 : 0.0;
      out["terz"] = (std::abs(fresnelQS[0]) + std::abs(fresnelQS[1]) > 1.0e-10) ? 
                    fresnelQS[0] * termRZ1 + fresnelQS[1] * termRZ2 : 0.0;
      out["tmrz"] = (std::abs(fresnelQS[2]) + std::abs(fresnelQS[3]) > 1.0e-10) ? 
                    fresnelQS[2] * termRZ1 + fresnelQS[3] * termRZ2 : 0.0;
    }
  }
  return out;
}

// Evaluate quasistatic Green's function components for **inter-layer** interactions.
std::map<std::string, dcmplx> GreenFLayered3D::evalGreenQSInter(const rvec pos1, const rvec pos2, 
                                                                double rho, double Z, double R,
                                                                bool singCancFields)
{    
  std::map<std::string,dcmplx> out; // Prepare the output map
  std::vector<dcmplx> fresnelQS;    // Quasistatic Fresnel reflection coefficients

  double dir = 0.0; // Propagation direction of reflected waves
  if (this->GetLayerIndex(pos1) < this->GetLayerIndex(pos2)) {
    dir = -1.0;
  }
  else {
    dir = 1.0;
  } 
  if (std::abs(this->GetLayerIndex(pos1) - this->GetLayerIndex(pos2)) == 1) {
    std::vector<dcmplx> temp = FresnelCoeff(this->GetLayerIndex(pos2), this->GetLayerIndex(pos1), 
                                            dcmplx(10000000000.0, 0.0), "TE");
    fresnelQS.push_back(temp[1]); // Transmission coefficient for TE wave
    temp = FresnelCoeff(this->GetLayerIndex(pos2), this->GetLayerIndex(pos1), 
                        dcmplx(10000000000.0, 0.0), "TM");
    fresnelQS.push_back(temp[1]); // Transmission coefficient for TM wave
  }
  else {
    fresnelQS.push_back(dcmplx(0.0, 0.0)); // Transmission coefficient for TE wave
    fresnelQS.push_back(dcmplx(0.0, 0.0)); // Transmission coefficient for TM wave
  }

  // quasistatic Green's functions terms
  dcmplx termZ  =  dir * log(Z + R) / (4.0 * PI);
  dcmplx termZZ = -1.0 / (4.0 * PI * R);
  dcmplx termS  =  1.0 / (4.0 * PI * R);
  dcmplx termR  = -(R - Z) / (4.0 * PI * rho);
  dcmplx termRZ =  dir * (1.0 - Z / R) / (4.0 * PI * rho);

  // form depending on singularity cancellation for E, H fields calculation
  if (singCancFields) {
    out["z1"]  =  termZ;
    out["z2"]  = -termZ;
    out["zz"]  =  termZZ;
    out["s"]   =  termS;
    out["r"]   =  termR;
    out["rz1"] =  termRZ;
    out["rz2"] = -termRZ;
    out["teTrans"] = fresnelQS[0];
    out["tmTrans"] = fresnelQS[1];
  } else {      
    out["tez1"]  = (std::abs(fresnelQS[0]) > 1.0e-10) ?  fresnelQS[0] * termZ : 0.0;
    out["tmz1"]  = (std::abs(fresnelQS[1]) > 1.0e-10) ?  fresnelQS[1] * termZ : 0.0;
    out["tez2"]  = (std::abs(fresnelQS[0]) > 1.0e-10) ? -fresnelQS[0] * termZ : 0.0;
    out["tmz2"]  = (std::abs(fresnelQS[1]) > 1.0e-10) ? -fresnelQS[1] * termZ : 0.0;
    out["tezz"]  = (std::abs(fresnelQS[0]) > 1.0e-10) ?  fresnelQS[0] * termZZ : 0.0;
    out["tmzz"]  = (std::abs(fresnelQS[1]) > 1.0e-10) ?  fresnelQS[1] * termZZ : 0.0;
    out["tes"]   = (std::abs(fresnelQS[0]) > 1.0e-10) ?  fresnelQS[0] * termS : 0.0;
    out["tms"]   = (std::abs(fresnelQS[1]) > 1.0e-10) ?  fresnelQS[1] * termS : 0.0;
    out["ter"]   = (std::abs(fresnelQS[0]) > 1.0e-10) ?  fresnelQS[0] * termR : 0.0;
    out["tmr"]   = (std::abs(fresnelQS[1]) > 1.0e-10) ?  fresnelQS[1] * termR : 0.0;
    out["terz1"] = (std::abs(fresnelQS[0]) > 1.0e-10) ?  fresnelQS[0] * termRZ : 0.0;
    out["tmrz1"] = (std::abs(fresnelQS[1]) > 1.0e-10) ?  fresnelQS[1] * termRZ : 0.0;
    out["terz2"] = (std::abs(fresnelQS[0]) > 1.0e-10) ? -fresnelQS[0] * termRZ : 0.0;
    out["tmrz2"] = (std::abs(fresnelQS[1]) > 1.0e-10) ? -fresnelQS[1] * termRZ : 0.0;
  }
  return out;
}

/**
 * @brief Generate uniformly spaced sample points in a 1D interval.
 * @param minVal Lower bound of the interval.
 * @param maxVal Upper bound of the interval.
 * @param numPoints Number of points to generate (>=2 for valid spacing).
 * @return Vector of uniformly spaced values from @p minVal to @p maxVal (inclusive).
 * @details Computes evenly spaced samples using a fixed step size 
 *          ( (maxVal - minVal) / (numPoints - 1) ).
 */
std::vector<double> generateUniformPoints(double minVal, double maxVal, int numPoints) {
  std::vector<double> points(numPoints);
  double step = (maxVal - minVal) / (numPoints - 1); // Uniform step size

  for (int i = 0; i < numPoints; ++i) {
      points[i] = minVal + i * step;
  }
  return points;
}

// Test function to fill the Green's function tables and write results to files.
void GreenFLayered3D::testFillGreenFTable(const std::string& intraTBType, const std::string& intraPMType, 
                                          const std::string& interType, const int& numPoints) {
  double x = 0.0, y = 0.0, z = 0.0;   // Values for interpolation results
  intraMap();  // Fill the function map for the intra-layer interactions
  interMap();  // Fill the function map for the inter-layer interactions
  // Counters for the member numbering of the interpolation coefficients structs
  int counterIntraTB = 0, counterIntraPM = 0, counterInter = 0;

  // Vectors for the interpolation results of each scenario
  std::vector<dcmplx> intraTBResults, intraPlusResults, 
                      intraMinusResults, interResults;

  // Iterate over all tabulation grids
  for (size_t i = 0; i < tabulationGrids.size(); i++) {
    // Top and bottom layers (if applicable)
    if ((tabulationGrids[i].i1 == tabulationGrids[i].i2) &&
        (tabulationGrids[i].i1 == 0 || tabulationGrids[i].i1 == epsilon.size() - 1) && 
        (intraTBIndex.size() > 0)) {
      
      // Find min and max values in the r vector
      auto minMaxR = std::minmax_element(tabulationGrids[i].r.begin(), tabulationGrids[i].r.end());
      double minValR = *minMaxR.first;
      double maxValR = *minMaxR.second;

      // Find min and max values in the z1 vector
      auto minMaxZ = std::minmax_element(tabulationGrids[i].z1.begin(), tabulationGrids[i].z1.end());
      double minValZ = *minMaxZ.first;
      double maxValZ = *minMaxZ.second;

      // Generate uniform points for r and z1
      std::vector<double> uniformR = generateUniformPoints(minValR, maxValR, numPoints);
      std::vector<double> uniformZ = generateUniformPoints(minValZ, maxValZ, numPoints);

      // Generate random r and z values and evaluate the Green's function for these
      for (int j = 0; j < numPoints; j++) {
        x = uniformR[j]; y = uniformZ[j];
        dcmplx res(this->evaluateSinglePoint(
                          x, y, coeffSelectorIntra[intraTBType](coeffIntraTB[counterIntraTB]), "real"), 
                   this->evaluateSinglePoint(
                          x, y, coeffSelectorIntra[intraTBType](coeffIntraTB[counterIntraTB]), "imag"));
        intraTBResults.push_back(res);
      }
      std::string fileName = "intraTBResults" + std::to_string(tabulationGrids[i].i1) + ".txt";
      this->writeVectorToFile(uniformR, uniformZ, intraTBResults, fileName);
      counterIntraTB++; intraTBResults.clear();
    }

    if ((tabulationGrids[i].i1 == tabulationGrids[i].i2) && 
        (tabulationGrids[i].i1 != 0 && tabulationGrids[i].i1 != epsilon.size() - 1) && 
        (intraPMIndex.size() > 0)) {
      
      // Find min and max values in the r vector
      auto minMaxR = std::minmax_element(tabulationGrids[i].r.begin(), tabulationGrids[i].r.end());
      double minValR = *minMaxR.first;
      double maxValR = *minMaxR.second;

      // Find min and max values in the z1 vector
      auto minMaxZ1 = std::minmax_element(tabulationGrids[i].z1.begin(), tabulationGrids[i].z1.end());
      double minValZ1 = *minMaxZ1.first;
      double maxValZ1 = *minMaxZ1.second;

      // Find min and max values in the z2 vector
      auto minMaxZ2 = std::minmax_element(tabulationGrids[i].z2.begin(), tabulationGrids[i].z2.end());
      double minValZ2 = *minMaxZ2.first;
      double maxValZ2 = *minMaxZ2.second;

      // Generate uniform points for r and z1
      std::vector<double> uniformR  = generateUniformPoints(minValR,  maxValR,  numPoints);
      std::vector<double> uniformZ1 = generateUniformPoints(minValZ1, maxValZ1, numPoints);
      std::vector<double> uniformZ2 = generateUniformPoints(minValZ2, maxValZ2, numPoints);

      // Generate random r and z values and evaluate the Green's function for these
      for (int j = 0; j < numPoints; j++) {
        x = uniformR[j]; y = uniformZ1[j]; z = uniformZ2[j];
        dcmplx res(this->evaluateSinglePoint(
                          x, y, coeffSelectorIntra[intraPMType](coeffIntraMinus[counterIntraPM]), "real"), 
                   this->evaluateSinglePoint(
                          x, y, coeffSelectorIntra[intraPMType](coeffIntraMinus[counterIntraPM]), "imag"));
        intraMinusResults.push_back(res);
        dcmplx res2(this->evaluateSinglePoint(
                          x, z, coeffSelectorIntra[intraPMType](coeffIntraPlus[counterIntraPM]), "real"), 
                    this->evaluateSinglePoint(
                          x, z, coeffSelectorIntra[intraPMType](coeffIntraPlus[counterIntraPM]), "imag"));
        intraPlusResults.push_back(res2);
      }
      std::string fileName1 = "intraMinusResults" + std::to_string(tabulationGrids[i].i1) + ".txt";
      std::string fileName2 = "intraPlusResults"  + std::to_string(tabulationGrids[i].i1) + ".txt";
      this->writeVectorToFile(uniformR, uniformZ1, intraMinusResults, fileName1);
      this->writeVectorToFile(uniformR, uniformZ2, intraPlusResults,  fileName2);
      counterIntraPM++; intraMinusResults.clear(); intraPlusResults.clear();
    }

    if (tabulationGrids[i].i1 != tabulationGrids[i].i2 && interIndeces.size() > 0) {
      
      // Find min and max values in the r vector
      auto minMaxR = std::minmax_element(tabulationGrids[i].r.begin(), tabulationGrids[i].r.end());
      double minValR = *minMaxR.first;
      double maxValR = *minMaxR.second;

      // Find min and max values in the z1 vector
      auto minMaxZ1 = std::minmax_element(tabulationGrids[i].z1.begin(), tabulationGrids[i].z1.end());
      double minValZ1 = *minMaxZ1.first;
      double maxValZ1 = *minMaxZ1.second;

      // Find min and max values in the z2 vector
      auto minMaxZ2 = std::minmax_element(tabulationGrids[i].z2.begin(), tabulationGrids[i].z2.end());
      double minValZ2 = *minMaxZ2.first;
      double maxValZ2 = *minMaxZ2.second;

      // Generate uniform points for r and z1
      std::vector<double> uniformR  = generateUniformPoints(minValR,  maxValR,  numPoints);
      std::vector<double> uniformZ1 = generateUniformPoints(minValZ1, maxValZ1, numPoints);
      std::vector<double> uniformZ2 = generateUniformPoints(minValZ2, maxValZ2, numPoints);

      // Generate random r and z values and evaluate the Green's function for these
      for (int j = 0; j < numPoints; j++) {
        x = uniformR[j]  * coeffInter[counterInter].scales[0]; 
        y = uniformZ1[j] * coeffInter[counterInter].scales[1];
        z = uniformZ2[j] * coeffInter[counterInter].scales[2];
        dcmplx res(this->evaluateSinglePoint(
                          x, y, z, coeffSelectorInter[interType](coeffInter[counterInter]), "real"), 
                   this->evaluateSinglePoint(
                          x, y, z, coeffSelectorInter[interType](coeffInter[counterInter]), "imag"));
        interResults.push_back(res);
      }
      std::string fileName = "interResults_obs" + std::to_string(tabulationGrids[i].i1) + "_src" + 
                                                  std::to_string(tabulationGrids[i].i2) + ".txt";
      this->writeVectorToFile(uniformR, uniformZ1, uniformZ2, interResults, fileName);
      counterInter++; interResults.clear();
    }
  }
}

// Write vector along with 2D coordinates to .txt file for debugging purposes.
void GreenFLayered3D::writeVectorToFile(const std::vector<double>& x, 
                                        const std::vector<double>& y, 
                                        const std::vector<dcmplx>& inVec, 
                                        const std::string& filename) {
  std::ofstream file(filename);
  if (!file) {
      std::cerr << "Error opening file: " << filename << std::endl;
      return;
  }

  for (size_t i = 0; i < inVec.size(); i++) {
      file << x[i] << " " << y[i] << " " 
           << real(inVec[i]) << " " << imag(inVec[i]) << "\n";
  }
  file.close();
}

// Write vector along with 2D coordinates to .txt file for debugging purposes.
void GreenFLayered3D::writeVectorToFile(const std::vector<rvec>& rho, 
                                        const std::vector<double>& inVec, 
                                        const std::string& filename) {
  std::ofstream file(filename);
  if (!file) {
      std::cerr << "Error opening file: " << filename << std::endl;
      return;
  }
  file << "A = [\n" << std::fixed << std::setprecision(6);
  for (size_t i = 0; i < inVec.size(); i++) {
      file << rho[i][0] << " " << rho[i][1] << " " << inVec[i] << "; ";
  }
  file << "];\n";
  file.close();
}

// Write vector along with 3D coordinates to .txt file for debugging purposes.
void GreenFLayered3D::writeVectorToFile(const std::vector<double>& x, 
                                        const std::vector<double>& y, 
                                        const std::vector<double>& z, 
                                        const std::vector<dcmplx>& inVec, 
                                        const std::string& filename) {
  std::ofstream file(filename);
  if (!file) {
      std::cerr << "Error opening file: " << filename << std::endl;
      return;
  }

  for (size_t i = 0; i < inVec.size(); i++) {
      file << x[i] << " " << y[i] << " " << z[i] << " " 
           << real(inVec[i]) << " " << imag(inVec[i]) << "\n";
  }
  file.close();
}

// Test the assembly process.
void GreenFLayered3D::testSameTriDK(std::vector<rvec> verts1, std::vector<rvec> verts2, 
                                    std::string gType, std::string vertPrint) {
  // Define observation and source triangle (same triangle)
  rvec barycenter1 = (verts1[0] + verts1[1] + verts1[2]) / 3.0;
  rvec barycenter2 = (verts2[0] + verts2[1] + verts2[2]) / 3.0;

  // Evaluate quadrature points
  int Ndunavant;         Ndunavant = dunavant::N17;
  double(*xdunavant)[3]; xdunavant = dunavant::x17;
  std::vector<rvec> svec(Ndunavant), spvec(Ndunavant);
  for (int i = 0; i < Ndunavant; ++i) {
    svec[i]  = verts1[0] * xdunavant[i][0] + verts1[1] * xdunavant[i][1] + verts1[2] * xdunavant[i][2];
    spvec[i] = verts2[0] * xdunavant[i][0] + verts2[1] * xdunavant[i][1] + verts2[2] * xdunavant[i][2];
  }

  // Set up the interaction type and coefficient index
  std::string interactionType = "unknown";
  int coeffIndex = -1;
  // Define the needed parameters for intra-layer and inter-layer interactions:
  // Source and observation layers are the same and are the top or bottom layer
  if ((GetLayerIndex(barycenter1) == GetLayerIndex(barycenter2)) &&
      (GetLayerIndex(barycenter1) == 0 || GetLayerIndex(barycenter1) == epsilon.size() - 1)) {
    interactionType = "intraTB";
    auto it = std::find(intraTBIndex.begin(), intraTBIndex.end(), GetLayerIndex(barycenter1));
    if (it == intraTBIndex.end()) throw std::runtime_error("Invalid coeffIndex in testSameTriDK");
    coeffIndex = std::distance(intraTBIndex.begin(), it);
  }
  // Source and observation layers are the same but they are not the first or last layer
  else if ((GetLayerIndex(barycenter1) == GetLayerIndex(barycenter2)) &&
          (GetLayerIndex(barycenter1) != 0 && GetLayerIndex(barycenter1) != epsilon.size() - 1)) {
    interactionType = "intraIn";
    auto it = std::find(intraPMIndex.begin(), intraPMIndex.end(), GetLayerIndex(barycenter1));
    if (it == intraPMIndex.end()) throw std::runtime_error("Invalid coeffIndex in testSameTriDK");
    coeffIndex = std::distance(intraPMIndex.begin(), it);
  }
  // Source and observation layers are different
  else if (GetLayerIndex(barycenter1) != GetLayerIndex(barycenter2)){
    interactionType = "inter";
    std::complex<int> target(GetLayerIndex(barycenter1), GetLayerIndex(barycenter2));
    auto it = std::find(interIndeces.begin(), interIndeces.end(), target);
    if (it == interIndeces.end()) throw std::runtime_error("Invalid coeffIndex in testSameTriDK");
    coeffIndex = std::distance(interIndeces.begin(), it);
  }
  else {
    std::cerr << "Something went wrong with the layers' indices during assembly evaluation" << std::endl;
  }

  if (vertPrint == "yes") {
    // Prepare output  for observation points
    std::ofstream outFile1("vertices1.txt");
    if (!outFile1) {
      std::cerr << "Error opening file: vertices1.txt" << std::endl;
      return;
    }
    outFile1 << "# x y z\n";

    for (int i = 0; i < Ndunavant; ++i) {
      outFile1 << svec[i][0] << " " << svec[i][1] << " " << svec[i][2] << ";\n";
    }
    outFile1.close();

    // Prepare output for source points
    std::ofstream outFile2("vertices2.txt");
    if (!outFile2) {
      std::cerr << "Error opening file: vertices2.txt" << std::endl;
      return;
    }
    outFile2 << "# x y z\n";

    for (int i = 0; i < Ndunavant; ++i) {
      outFile2 << spvec[i][0] << " " << spvec[i][1] << " " << spvec[i][2] << ";\n";
    }
    outFile2.close();    
  }
  
  // Prepare output
  std::cout << "Testing the calculation of smooth & quasistatic Green's function terms!" << std::endl;
  std::string filename = "test_SameTriDK_" + gType + ".txt";
  std::ifstream fileCheck(filename);
  if (!fileCheck) {
    std::ofstream outFile(filename);
    if (!outFile) {
      std::cerr << "Error opening file: " << filename << std::endl;
      return;
    }
    outFile << "# i j gTab_" << gType << " gQS_" << gType << " gTot_" << gType << "\n";

    for (int i = 0; i < Ndunavant; ++i) {
      for (int j = 0; j < Ndunavant; ++j) {
        auto coordSmooth = 
              this->layeredUtils->cartesianToLayeredTab(svec[i], spvec[j], interactionType);
        auto coordQS = 
              this->layeredUtils->cartesianToLayeredQS(svec[i], spvec[j], interactionType);
        auto gSmooth = (interactionType == "inter")
            ? this->evalGreenSmoothInter(coordSmooth.r, coordSmooth.Z1, coordSmooth.Z2, coeffIndex)
            : this->evalGreenSmoothIntra(interactionType, coordSmooth.r, coordSmooth.Z1, coordSmooth.Z2, coeffIndex);
        auto gQS = (interactionType == "inter")
            ? this->evalGreenQSInter(svec[i], spvec[j], coordQS.r, coordQS.Z1, coordQS.R1)
            : this->evalGreenQSIntra(svec[i], interactionType, coordQS.r, coordQS.Z1, 
                                      coordQS.Z2, coordQS.R1, coordQS.R2);

        outFile << i << " " << j << " " << gSmooth[gType] << " " 
                << gQS[gType] << " " << gSmooth[gType] + gQS[gType] << "\n";
      }
    }
    outFile.close();
  }
}