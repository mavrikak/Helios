/**
 * @file SommerfeldIntegrator.cpp
 * @brief Implementation of SommerfeldIntegrator: contour construction,
 *        integrand assembly, QUADPACK integration, and interpolation packing.
 */

#include "SommerfeldIntegrator.h"
#include <algorithm>
#include <numeric>
#include <unordered_map>

/* ----------------------------
 * File-scope state for QUADPACK callbacks
 * (QUADPACK expects function pointers; we memoize per-abscissa results
 * and use globals to select which member (TE/TM/derivatives) to return.)
 * ----------------------------
 * File-static utility variables for the Gauss-Kronrod integration, 
 * since no lambda function can be called as a pointer if it 
 * captures references in c++11:
 * ----------------------------
 */
static SommerfeldIntegrator* g_integ = nullptr; //!< Pointer to the SommerfeldIntegrator object (this)
static std::vector<dcmplx> g_rq_te;             //!< Quasistatic reflection coefficient (TE)
static std::vector<dcmplx> g_rq_tm;             //!< Quasistatic reflection coefficient (TM)
static dcmplx g_tq_te;                          //!< Quasistatic transmission coefficient (TE)
static dcmplx g_tq_tm;                          //!< Quasistatic transmission coefficient (TM)
static Data2D g_Data2D;                         //!< 2D data structure
static Data3D g_Data3D;                         //!< 3D data structure
static IntegrandIntra g_intZeroKrho;            //!< Integrand for intra-layer integrals at @f$k_\rho = 0@f$
static IntegrandInter g_intZeroKrhoInter;       //!< Integrand for inter-layer integrals at @f$k_\rho = 0@f$
static double g_k_max = 0.0;                    //!< Value for integration path definition
static std::string g_path = "none";             //!< Path of integration parameter (semi, real, imag)
static size_t g_currentI = 0;                   //!< Row of IntegrandIntra/IntegrandInter integration result
static size_t g_currentJ = 0;                   //!< Column of IntegrandIntra/IntegrandInter integration result
static size_t g_currentK = 0;                   //!< Depth of IntegrandIntra/IntegrandInter integration result
static std::string g_op  = "";                  //!< Operation type (minus, plus, "")

//! Enumeration of intra-layer integrand components
enum class IntraMember {
    TE, TM, TEZ, TMZ,
    TEZZ, TMZZ, TES, TMS,
    TER, TMR, TERZ, TMRZ
};

//! Enumeration of inter-layer integrand components
enum class InterMember {
    TE, TM, TEZ1, TMZ1, TEZ2, TMZ2,
    TEZZ, TMZZ, TES, TMS, TER, TMR, 
    TERZ1, TMRZ1, TERZ2, TMRZ2
};

static IntraMember g_currentMember = IntraMember::TE;      //!< Current intra-layer member (random initialization)
static InterMember g_currentMemberInter = InterMember::TE; //!< Current inter-layer member (random initialization)

//! Global memo cache: x -> IntegrandIntra
static std::map<double, IntegrandIntra> g_memo;

//! Global memo cache: x -> IntegrandInter
static std::unordered_map<double, IntegrandInter> g_memoInter;

//! Intra or inter global status
std::string g_status = "";

// ----------------------------
// QUADPACK callback selectors
// ----------------------------
/**
 * @brief Retrieves the real part of a selected member of the computed integrand.
 *
 * This function selects the correct member from either the "inter" or "intra" integrand,
 * computing it if necessary and caching the result for future lookups.
 *
 * @param x Pointer to the independent variable for which the integrand is evaluated.
 * @return The real part of the selected member of the computed integrand.
 *
 * @throws std::invalid_argument if g_status is neither "inter" nor "intra".
 */
extern "C" double pickRealMember(double* x)
{   
    dcmplx c(0.0, 0.0);
    if (g_status == "inter") {
        // Check if we already have IntegrandInter computed for x
        auto it = g_memoInter.find(*x);
        if (it == g_memoInter.end())
        {
            // Not found in cache -> Evaluate the entire integrand at x
            IntegrandInter val = g_integ->integrationFInter(
                g_tq_te, g_tq_tm, g_Data3D, *x, g_k_max, g_intZeroKrhoInter, g_path
            );

            // Store in cache
            it = g_memoInter.emplace(*x, std::move(val)).first;
        }

        // Now "it->second" is the cached IntegrandInter
        const IntegrandInter &val = it->second;

        // Pick the correct (member, i, j, k)
        switch (g_currentMemberInter) {
            case InterMember::TE:    c = val.te[g_currentI][g_currentJ][g_currentK];    break;
            case InterMember::TM:    c = val.tm[g_currentI][g_currentJ][g_currentK];    break;
            case InterMember::TEZ1:  c = val.tez1[g_currentI][g_currentJ][g_currentK];  break;
            case InterMember::TMZ1:  c = val.tmz1[g_currentI][g_currentJ][g_currentK];  break;
            case InterMember::TEZ2:  c = val.tez2[g_currentI][g_currentJ][g_currentK];  break;
            case InterMember::TMZ2:  c = val.tmz2[g_currentI][g_currentJ][g_currentK];  break;
            case InterMember::TEZZ:  c = val.tezz[g_currentI][g_currentJ][g_currentK];  break;
            case InterMember::TMZZ:  c = val.tmzz[g_currentI][g_currentJ][g_currentK];  break;
            case InterMember::TES:   c = val.tes[g_currentI][g_currentJ][g_currentK];   break;
            case InterMember::TMS:   c = val.tms[g_currentI][g_currentJ][g_currentK];   break;
            case InterMember::TER:   c = val.ter[g_currentI][g_currentJ][g_currentK];   break;
            case InterMember::TMR:   c = val.tmr[g_currentI][g_currentJ][g_currentK];   break;
            case InterMember::TERZ1: c = val.terz1[g_currentI][g_currentJ][g_currentK]; break;
            case InterMember::TMRZ1: c = val.tmrz1[g_currentI][g_currentJ][g_currentK]; break;
            case InterMember::TERZ2: c = val.terz2[g_currentI][g_currentJ][g_currentK]; break;
            case InterMember::TMRZ2: c = val.tmrz2[g_currentI][g_currentJ][g_currentK]; break;
        }
    } else if (g_status == "intra"){
        // Check if we already have IntegrandIntra computed for x
        auto it = g_memo.find(*x);
        if (it == g_memo.end())
        {
            // Not found in cache -> Evaluate the entire integrand at x
            IntegrandIntra val = g_integ->integrationFIntra(
                g_rq_te, g_rq_tm, g_Data2D, *x, g_k_max, g_intZeroKrho, g_path, g_op
            );

            // Store in the map
            g_memo[*x] = val;
            it = g_memo.find(*x); // "it" now points to the newly inserted element
        }

        // Now "it->second" is the cached IntegrandIntra
        const IntegrandIntra &val = it->second;

        // Pick the correct (member, i, j)
        switch (g_currentMember) {
            case IntraMember::TE:   c = val.te[g_currentI][g_currentJ];   break;
            case IntraMember::TM:   c = val.tm[g_currentI][g_currentJ];   break;
            case IntraMember::TEZ:  c = val.tez[g_currentI][g_currentJ];  break;
            case IntraMember::TMZ:  c = val.tmz[g_currentI][g_currentJ];  break;
            case IntraMember::TEZZ: c = val.tezz[g_currentI][g_currentJ]; break;
            case IntraMember::TMZZ: c = val.tmzz[g_currentI][g_currentJ]; break;
            case IntraMember::TES:  c = val.tes[g_currentI][g_currentJ];  break;
            case IntraMember::TMS:  c = val.tms[g_currentI][g_currentJ];  break;
            case IntraMember::TER:  c = val.ter[g_currentI][g_currentJ];  break;
            case IntraMember::TMR:  c = val.tmr[g_currentI][g_currentJ];  break;
            case IntraMember::TERZ: c = val.terz[g_currentI][g_currentJ]; break;
            case IntraMember::TMRZ: c = val.tmrz[g_currentI][g_currentJ]; break;
        }
    } else {
        throw std::invalid_argument("Invalid status. Use 'intra' or 'inter'.");
    }
    return std::real(c);
}

/**
 * @brief Retrieves the imaginary part of a selected member of the computed integrand.
 *
 * This function selects the correct member from either the "inter" or "intra" integrand,
 * computing it if necessary and caching the result for future lookups.
 *
 * @param x Pointer to the independent variable for which the integrand is evaluated.
 * @return The imaginary part of the selected member of the computed integrand.
 *
 * @throws std::invalid_argument if g_status is neither "inter" nor "intra".
 */
extern "C" double pickImagMember(double* x)
{
    dcmplx c(0.0, 0.0);
    if (g_status == "inter") {
        // Check if we already have IntegrandInter computed for x
        auto it = g_memoInter.find(*x);
        if (it == g_memoInter.end())
        {
            // Not found in cache -> Evaluate the entire integrand at x
            IntegrandInter val = g_integ->integrationFInter(
                g_tq_te, g_tq_tm, g_Data3D, *x, g_k_max, g_intZeroKrhoInter, g_path
            );

            // Store in cache
            it = g_memoInter.emplace(*x, std::move(val)).first;
        }

        // Now "it->second" is the cached IntegrandIntra
        const IntegrandInter &val = it->second;

        // Pick the correct (member, i, j, k)
        switch (g_currentMemberInter) {
            case InterMember::TE:    c = val.te[g_currentI][g_currentJ][g_currentK];    break;
            case InterMember::TM:    c = val.tm[g_currentI][g_currentJ][g_currentK];    break;
            case InterMember::TEZ1:  c = val.tez1[g_currentI][g_currentJ][g_currentK];  break;
            case InterMember::TMZ1:  c = val.tmz1[g_currentI][g_currentJ][g_currentK];  break;
            case InterMember::TEZ2:  c = val.tez2[g_currentI][g_currentJ][g_currentK];  break;
            case InterMember::TMZ2:  c = val.tmz2[g_currentI][g_currentJ][g_currentK];  break;
            case InterMember::TEZZ:  c = val.tezz[g_currentI][g_currentJ][g_currentK];  break;
            case InterMember::TMZZ:  c = val.tmzz[g_currentI][g_currentJ][g_currentK];  break;
            case InterMember::TES:   c = val.tes[g_currentI][g_currentJ][g_currentK];   break;
            case InterMember::TMS:   c = val.tms[g_currentI][g_currentJ][g_currentK];   break;
            case InterMember::TER:   c = val.ter[g_currentI][g_currentJ][g_currentK];   break;
            case InterMember::TMR:   c = val.tmr[g_currentI][g_currentJ][g_currentK];   break;
            case InterMember::TERZ1: c = val.terz1[g_currentI][g_currentJ][g_currentK]; break;
            case InterMember::TMRZ1: c = val.tmrz1[g_currentI][g_currentJ][g_currentK]; break;
            case InterMember::TERZ2: c = val.terz2[g_currentI][g_currentJ][g_currentK]; break;
            case InterMember::TMRZ2: c = val.tmrz2[g_currentI][g_currentJ][g_currentK]; break;
        }
    } else if (g_status == "intra"){
        // Check if we already have IntegrandIntra computed for x
        auto it = g_memo.find(*x);
        if (it == g_memo.end())
        {
            // Not found in cache -> Evaluate the entire integrand at x
            IntegrandIntra val = g_integ->integrationFIntra(
                g_rq_te, g_rq_tm, g_Data2D, *x, g_k_max, g_intZeroKrho, g_path, g_op
            );

            // Store in the map
            g_memo[*x] = val;
            it = g_memo.find(*x); // "it" now points to the newly inserted element
        }

        // Now "it->second" is the cached IntegrandIntra
        const IntegrandIntra &val = it->second;

        // Pick the correct (member, i, j)
        switch (g_currentMember) {
            case IntraMember::TE:   c = val.te[g_currentI][g_currentJ];   break;
            case IntraMember::TM:   c = val.tm[g_currentI][g_currentJ];   break;
            case IntraMember::TEZ:  c = val.tez[g_currentI][g_currentJ];  break;
            case IntraMember::TMZ:  c = val.tmz[g_currentI][g_currentJ];  break;
            case IntraMember::TEZZ: c = val.tezz[g_currentI][g_currentJ]; break;
            case IntraMember::TMZZ: c = val.tmzz[g_currentI][g_currentJ]; break;
            case IntraMember::TES:  c = val.tes[g_currentI][g_currentJ];  break;
            case IntraMember::TMS:  c = val.tms[g_currentI][g_currentJ];  break;
            case IntraMember::TER:  c = val.ter[g_currentI][g_currentJ];  break;
            case IntraMember::TMR:  c = val.tmr[g_currentI][g_currentJ];  break;
            case IntraMember::TERZ: c = val.terz[g_currentI][g_currentJ]; break;
            case IntraMember::TMRZ: c = val.tmrz[g_currentI][g_currentJ]; break;
        }
    } else {
        throw std::invalid_argument("Invalid status. Use 'intra' or 'inter'.");
    }
    return std::imag(c);
}

// ----------------------------
// External Fortran routines for Bessel/Hankel functions, 
// and Gauss-Kronrod numerical integration.
// ----------------------------
extern "C" {
    /**
     * @brief Compute complex Bessel functions of the first kind @f$ J_\nu(z) @f$.
     *
     * Fortran wrapper to the Amos library routine `ZBESJ`.
     *
     * @param zr Real part of the complex argument z.
     * @param zi Imaginary part of the complex argument z.
     * @param fnu Order of the Bessel function.
     * @param kode Scaling option (1 = unscaled, 2 = exponentially scaled).
     * @param n Number of terms to compute (number of sequential orders).
     * @param cyr Output array for the real parts of @f$ J_\nu(z) @f$.
     * @param cyi Output array for the imaginary parts of @f$ J_\nu(z) @f$.
     * @param nz Number of underflows set to zero.
     * @param ierr Error flag (0 if successful).
     */
    void zbesj_(const double* zr, const double* zi,
              const double* fnu, const int* kode, const int* n,
              double* cyr, double* cyi, int* nz, int* ierr);
    
    /**
     * @brief Compute complex Hankel functions of the first kind or second kind @f$ H_\nu^{(m)}(z) @f$.
     *
     * Fortran wrapper to the Amos library routine `ZBESH`.
     *
     * @param zr Real part of the complex argument z.
     * @param zi Imaginary part of the complex argument z.
     * @param fnu Order of the Hankel function.
     * @param kode Scaling option (1 = unscaled, 2 = exponentially scaled).
     * @param m Type selector (1 = first kind, 2 = second kind).
     * @param n Number of sequential orders to compute.
     * @param cyr Output array for the real parts of @f$ H_\nu^{(m)}(z) @f$.
     * @param cyi Output array for the imaginary parts of @f$ H_\nu^{(m)}(z) @f$.
     * @param nz Number of underflows set to zero.
     * @param ierr Error flag (0 if successful).
     */
    void zbesh_(const double* zr, const double* zi,
                const double* fnu, const int* kode, const int* m, const int* n,
                double* cyr, double* cyi, int* nz, int* ierr);
    
    /**
     * @brief Adaptive Gauss–Kronrod quadrature (QUADPACK driver routine DQAGE).
     *
     * Evaluates a real integral using adaptive Gauss–Kronrod rules over a finite interval.
     *
     * @param f Pointer to the function to integrate.
     * @param a Lower integration bound.
     * @param b Upper integration bound.
     * @param epsabs Absolute accuracy requested.
     * @param epsrel Relative accuracy requested.
     * @param key Rule selection key (1–6 corresponding to different Gauss–Kronrod pairs).
     * @param limit Maximum number of subintervals allowed.
     * @param result Output integral result.
     * @param abserr Estimated absolute error.
     * @param neval Number of function evaluations performed.
     * @param ier Error flag (0 if successful).
     * @param alist Workspace array of left subinterval bounds.
     * @param blist Workspace array of right subinterval bounds.
     * @param rlist Workspace array of subinterval results.
     * @param elist Workspace array of subinterval errors.
     * @param iord Workspace array for ordering subintervals.
     * @param last Number of subintervals actually used.
     */
    void dqage_(double (*f)(double *), double *a, double *b,
                double *epsabs, double *epsrel, int *key, int *limit,
                double *result, double *abserr, int *neval, int *ier,
                double *alist, double *blist, double *rlist, double *elist,
                int *iord, int *last);
}

// ----------------------------
// Special functions / integrator wrapper
// ----------------------------
/**
 * @brief Computes the Bessel function @f$ J_n(z) @f$.
 * @param order Order of the Bessel function.
 * @param arg Argument of the function.
 * @return Computed Bessel function value.
 */
dcmplx besselJ(int order, dcmplx arg) {
    // Large-argument shortcut
    if (std::abs(arg) >= 32768.0) {
        return sqrt(2.0 / (M_PI * arg)) * cos(arg - order * M_PI / 2.0 - M_PI / 4.0);
    }

    const double zr = std::real(arg), zi = std::imag(arg);
    const double fnu = static_cast<double>(order);
    const int n = 1, kode = 1; // unscaled
    double cyr = 0.0, cyi = 0.0;
    int nz = 0, ierr = 0;

    zbesj_(&zr, &zi, &fnu, &kode, &n, &cyr, &cyi, &nz, &ierr);
    if (ierr != 0) {
        std::cerr << "AMOS ZBESJ error: ierr=" << ierr
                  << " for J_" << order << "(" << arg << ")\n";
    }
    return dcmplx(cyr, cyi);
}

/**
 * @brief Computes the Hankel function @f$ H_\nu^{(m)}(z) @f$.
 * @param order Order of the Hankel function.
 * @param kind Kind of Hankel function (1 or 2).
 * @param arg Argument of the function.
 * @return Computed Hankel function value.
 */
dcmplx hankelH(int order, int kind, dcmplx arg) {
    // Large-argument shortcut
    if (std::abs(arg) >= 32768.0) {
        if (kind == 1) {
            return sqrt(2.0 / (M_PI * arg)) * exp( I * (arg - order * M_PI / 2.0 - M_PI / 4.0));
        } else { // kind == 2
            return sqrt(2.0 / (M_PI * arg)) * exp(-I * (arg - order * M_PI / 2.0 - M_PI / 4.0));
        }
    }

    const double zr = std::real(arg), zi = std::imag(arg);
    const double fnu = static_cast<double>(order);
    const int n = 1, kode = 1;
    const int m = (kind == 2 ? 2 : 1); // 1 -> H^(1), 2 -> H^(2)
    double cyr = 0.0, cyi = 0.0;
    int nz = 0, ierr = 0;

    zbesh_(&zr, &zi, &fnu, &kode, &m, &n, &cyr, &cyi, &nz, &ierr);
    if (ierr != 0) {
        std::cerr << "AMOS ZBESH error: ierr=" << ierr
                  << " for H_" << (m==1?"(1)":"(2)") << order << "(" << arg << ")\n";
    }
    return dcmplx(cyr, cyi);
}

/**
 * @brief Performs Gauss-Kronrod integration.
 * @param integrand Function pointer to the integrand function.
 * @param a Lower bound of integration.
 * @param b Upper bound of integration.
 * @param key Selection of Gauss-Kronrod pair (5 -> 25/51).
 * @param limit Maximum subdivisions.
 * @param epsabs Absolute error tolerance.
 * @param epsrel Relative error tolerance.
 * @return Computed integral value as a complex number.
 */
dcmplx gaussKronrod(double (*integrand)(double*), double a, double b,
                    int key = 5, int limit = 100,
                    double epsabs = 1e-6, double epsrel = 1e-6)
{
    // QUADPACK expects key in [1..6] -> 15,21,31,41,51,61 (51 <-> 25/51 pair)
    int KEY   = std::max(1, std::min(6, key));
    int LIMIT = std::max(1, limit);

    // Workspaces
    std::vector<double> ALIST(LIMIT), BLIST(LIMIT), RLIST(LIMIT), ELIST(LIMIT);
    std::vector<int>    IORD(LIMIT);

    // Outputs
    double result = 0.0, abserr = 0.0;
    int neval = 0, last = 0, ier = 0;

    dqage_(integrand, &a, &b, &epsabs, &epsrel, &KEY, &LIMIT,
           &result, &abserr, &neval, &ier,
           ALIST.data(), BLIST.data(), RLIST.data(), ELIST.data(),
           IORD.data(), &last);

    if (ier != 0) {
        std::cerr << "[QUADPACK DQAG] ier=" << ier
                  << "  result=" << result
                  << "  abserr=" << abserr
                  << "  neval=" << neval
                  << "  last="  << last << "\n";
    }

    // Return result with error in imag part
    return dcmplx(result, abserr);
}

// ----------------------------
// Interpolation helpers (packing)
// ----------------------------
/**
 * @brief Flatten 2D complex grid to column-major array with axis metadata.
 * @param x_grid A vector containing the x-coordinates of the grid.
 * @param y_grid A vector containing the y-coordinates of the grid.
 * @param green A 2D vector containing the complex values of the function on the grid.
 * @return A tuple containing:
 *         - F_flat: green flattened in column-major order.
 *         - px: A vector of two integers representing px values for real and imaginary parts.
 *         - py: A vector of two integers representing py values for real and imaginary parts.
 *         - Xc: x_grid points.
 *         - Yc: y_grid points.
 */
std::tuple<std::vector<dcmplx>, std::vector<int>, std::vector<int>,
           std::vector<double>, std::vector<double>>
interpolateIntra(std::vector<double>& x_grid, std::vector<double>& y_grid,
                 std::vector<std::vector<dcmplx>>& green)
{
    const int mx = static_cast<int>(x_grid.size());
    const int my = static_cast<int>(y_grid.size());

    // Flatten F to Fortran column-major: first index (x) runs fastest
    std::vector<dcmplx> F_flat(mx * my);
    for (int j = 0; j < my; ++j)            // y index (slow)
      for (int i = 0; i < mx; ++i)          // x index (fast)
        F_flat[i + j*mx] = green[i][j];

    // Carry grid boundaries as complex (imag=0) to fit existing tuple type
    std::vector<double> Xc(mx), Yc(my);
    for (int i = 0; i < mx; ++i) Xc[i] = x_grid[i];
    for (int j = 0; j < my; ++j) Yc[j] = y_grid[j];

    // px/py vectors store counts (mx,my) for both real/imag branches
    std::vector<int> px{mx, mx}, py{my, my};

    return std::make_tuple(std::move(F_flat), std::move(px), std::move(py),
                           std::move(Xc), std::move(Yc));
}

/**
 * @brief Flatten 3D complex tensor grid (column-major order).
 * @param x A 3D vector representing the x-coordinates of the grid.
 * @param y A 3D vector representing the y-coordinates of the grid.
 * @param z A 3D vector representing the z-coordinates of the grid.
 * @param inScales A vector with the scaling factors for x, y, z coordinates.
 * @param green A 3D vector containing complex function values at grid points.
 * @return A tuple containing:
 *         - Xv: A vector of x-coordinates.
 *         - Yv: A vector of y-coordinates.
 *         - Zv: A vector of z-coordinates.
 *         - F_flat: green flattened in column-major order.
 */
std::tuple<std::vector<double>, std::vector<double>, std::vector<double>, std::vector<dcmplx>>
interpolateInter(std::vector<std::vector<std::vector<double>>>& x,
                 std::vector<std::vector<std::vector<double>>>& y,
                 std::vector<std::vector<std::vector<double>>>& z,
                 std::vector<double>& inScales,
                 std::vector<std::vector<std::vector<dcmplx>>>& green)
{
    if (x.empty() || x[0].empty() || x[0][0].empty())
        throw std::runtime_error("Error: Input arrays for interpolation cannot be empty.");

    const int mx = static_cast<int>(x.size());
    const int my = static_cast<int>(x[0].size());
    const int mz = static_cast<int>(x[0][0].size());

    // 1D axes (assumed rectilinear tensor grid)
    std::vector<double> Xv(mx), Yv(my), Zv(mz);
    for (int i = 0; i < mx; ++i) Xv[i] = inScales[0] * x[i][0][0];
    for (int j = 0; j < my; ++j) Yv[j] = inScales[1] * y[0][j][0];
    for (int k = 0; k < mz; ++k) Zv[k] = inScales[2] * z[0][0][k];

    // Column-major flattening: i (x) fastest, then j (y), then k (z)
    std::vector<dcmplx> F_flat(mx * my * mz);
    for (int k = 0; k < mz; ++k)
      for (int j = 0; j < my; ++j)
        for (int i = 0; i < mx; ++i)
          F_flat[i + j*mx + k*mx*my] = green[i][j][k];

    return std::make_tuple(std::move(Xv), std::move(Yv), std::move(Zv), std::move(F_flat));
}

// ------------------------------------------------
// Constructors & Destructor
// ------------------------------------------------
SommerfeldIntegrator::SommerfeldIntegrator(const std::vector<dcmplx> inLayeredK, int inMaterialIndex,
                                           const std::vector<dcmplx>& inEpsilonMedium,
                                           const std::vector<dcmplx>& inMuMedium,
                                           const std::vector<double>& inZValsInterfaces,
                                           const std::vector<double>& inThickness,
                                           double inSemiReal, double inSemiImag,
                                           double inRatio, double inCutoff)
    : k_L(inLayeredK), materialIndex(inMaterialIndex), epsilon(inEpsilonMedium), mu(inMuMedium),
      zValsInterfaces(inZValsInterfaces), thickness(inThickness), semiReal(inSemiReal),
      semiImag(inSemiImag), ratio(inRatio), cutoff(inCutoff) {
        layeredUtils = new LayeredMediaUtils(k_L, epsilon, mu, zValsInterfaces, thickness);
    }

SommerfeldIntegrator::SommerfeldIntegrator(const std::vector<dcmplx> inLayeredK, 
                                           int inMaterialIndex, int inSourceIndex,
                                           const std::vector<dcmplx>& inEpsilonMedium,
                                           const std::vector<dcmplx>& inMuMedium,
                                           const std::vector<double>& inZValsInterfaces,
                                           const std::vector<double>& inThickness,
                                           double inSemiReal, double inSemiImag,
                                           double inRatio, double inCutoff)
    : k_L(inLayeredK), materialIndex(inMaterialIndex), sourceIndex(inSourceIndex), 
      epsilon(inEpsilonMedium), mu(inMuMedium), zValsInterfaces(inZValsInterfaces), 
      thickness(inThickness), semiReal(inSemiReal), semiImag(inSemiImag), 
      ratio(inRatio), cutoff(inCutoff) {
        layeredUtils = new LayeredMediaUtils(k_L, epsilon, mu, zValsInterfaces, thickness);
    }

SommerfeldIntegrator::~SommerfeldIntegrator() {
    delete layeredUtils; // Free allocated memory
    layeredUtils = nullptr;
}

// ----------------------------
// Data selectors
// ----------------------------
/**
 * @details Pack 2D r/z data depending on path.
 *          - "semi": keep matrices; dense traversal by linear index.
 *          - "real": select (r,z) with  z >= ratio*r   (or z1+z2 >= ratio*r if provided).
 *          - "imag": select (r,z) with  z <  ratio*r   (or z1+z2 <  ratio*r if provided).
 */
Data2D SommerfeldIntegrator::selectData2D(const std::vector<std::vector<double>>& r,
                                          const std::vector<std::vector<double>>& z1,
                                          const std::vector<std::vector<double>>& z2,
                                          const std::string& path) {
    Data2D data;
    if (z2.empty()) {
        data.z2Vec.clear();
        data.z2Matrix.clear();
    }

    if (path == "semi") {
        data.rMatrix.resize(r.size(), std::vector<double>(r[0].size(), 0.0));
        data.z1Matrix.resize(z1.size(), std::vector<double>(z1[0].size(), 0.0));
        data.rMatrix = r; data.rVec.clear();
        data.z1Matrix = z1; data.z1Vec.clear();
        if (!z2.empty()) {
            data.z2Matrix.resize(z2.size(), std::vector<double>(z2[0].size(), 0.0));
            data.z2Matrix = z2; data.z2Vec.clear();
        }
        data.indVec.resize(r.size() * r[0].size());
        // Fill indices with 0, 1, ..., r.size() * r[0].size() - 1
        std::iota(data.indVec.begin(), data.indVec.end(), 0);
        data.indVecComplex.clear();
    } else if (path == "real") {
        data.rVec.clear(); data.z1Vec.clear(); data.z2Vec.clear();
        data.rMatrix.clear(); data.z1Matrix.clear(); data.z2Matrix.clear();
        data.indVec.clear(); data.indVecComplex.clear();
        for (size_t j = 0; j < r[0].size(); ++j) {
            for (size_t i = 0; i < r.size(); ++i) {
                if (z2.empty() && z1[i][j] >= this->ratio * r[i][j]) {
                    data.rVec.push_back(r[i][j]);
                    data.z1Vec.push_back(z1[i][j]);
                    data.indVecComplex.push_back(std::complex<int>(i, j));
                } else if (!z2.empty() && z1[i][j] + z2[i][j] >= this->ratio * r[i][j]) {
                    data.rVec.push_back(r[i][j]);
                    data.z1Vec.push_back(z1[i][j]);
                    data.z2Vec.push_back(z2[i][j]);
                    data.indVecComplex.push_back(std::complex<int>(i, j));
                }
            }
        }
    } else if (path == "imag") {
        data.rVec.clear(); data.z1Vec.clear(); data.z2Vec.clear();
        data.rMatrix.clear(); data.z1Matrix.clear(); data.z2Matrix.clear();
        data.indVec.clear(); data.indVecComplex.clear();
        for (size_t j = 0; j < r[0].size(); ++j) {
            for (size_t i = 0; i < r.size(); ++i) {
                if (z2.empty() && z1[i][j] < this->ratio * r[i][j]) {
                    data.rVec.push_back(r[i][j]);
                    data.z1Vec.push_back(z1[i][j]);
                    data.indVecComplex.push_back(std::complex<int>(i, j));
                } else if (!z2.empty() && z1[i][j] + z2[i][j] < this->ratio * r[i][j]) {
                    data.rVec.push_back(r[i][j]);
                    data.z1Vec.push_back(z1[i][j]);
                    data.z2Vec.push_back(z2[i][j]);
                    data.indVecComplex.push_back(std::complex<int>(i, j));
                }
            }
        }
    } else {
        throw std::invalid_argument("Invalid path type. Use 'semi', 'real', or 'imag'.");
    }

    return data;
}

/**
 * @details Pack 3D r/z1/z2 data depending on path (semi/real/imag).
 *          The real/imag split uses |z1 - z2| vs ratio*r.
 */
Data3D SommerfeldIntegrator::selectData3D(
            const std::vector<std::vector<std::vector<double>>>& r,
            const std::vector<std::vector<std::vector<double>>>& z1,
            const std::vector<std::vector<std::vector<double>>>& z2,
            const std::string& path) 
{
    Data3D data;
    if (path == "semi") {
        data.rMatrix.resize(r.size(), 
                            std::vector<std::vector<double>>(r[0].size(), 
                            std::vector<double>(r[0][0].size(), 0.0)));
        data.z1Matrix.resize(z1.size(), 
                             std::vector<std::vector<double>>(z1[0].size(), 
                             std::vector<double>(z1[0][0].size(), 0.0)));
        data.z2Matrix.resize(z2.size(), 
                             std::vector<std::vector<double>>(z2[0].size(), 
                             std::vector<double>(z2[0][0].size(), 0.0)));
        data.rMatrix = r;   data.rVec.clear();
        data.z1Matrix = z1; data.z1Vec.clear();
        data.z2Matrix = z2; data.z2Vec.clear();
        data.indVec.resize(r.size() * r[0].size() * r[0][0].size());
        // Fill indices with 0, 1, ..., r.size() * r[0].size() * r[0][0].size() - 1
        std::iota(data.indVec.begin(), data.indVec.end(), 0);
        data.indVecTriplet.clear();
    } else if (path == "real") {
        data.rVec.clear();    data.z1Vec.clear();    data.z2Vec.clear();
        data.rMatrix.clear(); data.z1Matrix.clear(); data.z2Matrix.clear();
        data.indVec.clear();  data.indVecTriplet.clear();
        for (size_t k = 0; k < r[0][0].size(); ++k) {
            for (size_t j = 0; j < r[0].size(); ++j) {
                for (size_t i = 0; i < r.size(); ++i) {
                    if (std::abs(z1[i][j][k] - z2[i][j][k]) >= this->ratio * r[i][j][k]) {
                        data.rVec.push_back(r[i][j][k]);
                        data.z1Vec.push_back(z1[i][j][k]);
                        data.z2Vec.push_back(z2[i][j][k]);
                        data.indVecTriplet.emplace_back(i, j, k);
                    }
                }
            }
        }
    } else if (path == "imag") {
        data.rVec.clear();    data.z1Vec.clear();    data.z2Vec.clear();
        data.rMatrix.clear(); data.z1Matrix.clear(); data.z2Matrix.clear();
        data.indVec.clear();  data.indVecTriplet.clear();
        for (size_t k = 0; k < r[0][0].size(); ++k) {
            for (size_t j = 0; j < r[0].size(); ++j) {
                for (size_t i = 0; i < r.size(); ++i) {
                    if (std::abs(z1[i][j][k] - z2[i][j][k]) < this->ratio * r[i][j][k]) {
                        data.rVec.push_back(r[i][j][k]);
                        data.z1Vec.push_back(z1[i][j][k]);
                        data.z2Vec.push_back(z2[i][j][k]);
                        data.indVecTriplet.emplace_back(i, j, k);
                    }
                }
            }
        }
    } else {
        throw std::invalid_argument("Invalid path type. Use 'semi', 'real', or 'imag'.");
    }

    return data;
}

// Shift z distances to nearest interface (sign via direction).
std::vector<std::vector<std::vector<double>>> 
SommerfeldIntegrator::interDist(std::vector<std::vector<std::vector<double>>> z, 
                                int index, double direction) {
    std::vector<std::vector<std::vector<double>>> result;
    size_t rows  = z.size();
    size_t cols  = z[0].size();
    size_t depth = z[0][0].size();

    if (direction == 1.0) {
        for (size_t i = 0; i < rows; i++) {
            for (size_t j = 0; j < cols; j++) {
                for (size_t k = 0; k < depth; k++) {
                    z[i][j][k] -= this->zValsInterfaces[index - 1];
                }
            }
        }
    } else {
        for (size_t i = 0; i < rows; i++) {
            for (size_t j = 0; j < cols; j++) {
                for (size_t k = 0; k < depth; k++) {
                    z[i][j][k] = this->zValsInterfaces[index] - z[i][j][k];
                }
            }
        }
    }
    return z;
}

// ----------------------------
// Secondary waves (inter-layer)
// ----------------------------
// Compose inter-layer secondary fields and their Z-derivatives for TE/TM.
std::vector<std::vector<std::vector<std::vector<dcmplx>>>> 
SommerfeldIntegrator::Secondary(dcmplx krho, 
                                std::vector<std::vector<std::vector<double>>> r,
                                std::vector<std::vector<std::vector<double>>> z1,
                                std::vector<std::vector<std::vector<double>>> z2,
                                const std::string& waves)
{   
    std::vector<std::vector<std::vector<dcmplx>>> f, fZ1, fZ2, fZZ;

    // Coefficients for TE/TM secondary waves
    blitz::Array<dcmplx, 2> waveSec = 
                        this->layeredUtils->Secondary(krho, this->materialIndex, this->sourceIndex, waves);
    
    // Wavenumber z-component for observation and source layers of inter-layer calculations
    dcmplx kz1 = this->k_L[this->materialIndex]; kz1 = sqrt(kz1 * kz1 - krho * krho);
    dcmplx kz2 = this->k_L[this->sourceIndex];   kz2 = sqrt(kz2 * kz2 - krho * krho);
    if (std::imag(kz1) < 0)  kz1 = -kz1; // Ensure positive imaginary part
    if (std::imag(kz2) < 0)  kz2 = -kz2; // Ensure positive imaginary part

    size_t rows  = r.size();
    size_t cols  = r[0].size();
    size_t depth = r[0][0].size();

    // Initialize the result structure components
    f.assign(  rows, std::vector<std::vector<dcmplx>>(cols, std::vector<dcmplx>(depth, dcmplx(0.0, 0.0) ) ) );
    fZ1.assign(rows, std::vector<std::vector<dcmplx>>(cols, std::vector<dcmplx>(depth, dcmplx(0.0, 0.0) ) ) );
    fZ2.assign(rows, std::vector<std::vector<dcmplx>>(cols, std::vector<dcmplx>(depth, dcmplx(0.0, 0.0) ) ) );
    fZZ.assign(rows, std::vector<std::vector<dcmplx>>(cols, std::vector<dcmplx>(depth, dcmplx(0.0, 0.0) ) ) );

    std::vector<double> dir = {1.0, -1.0}; // Direction of the waves

    for (int ii = 0; ii < waveSec.extent(0); ++ii) { // Rows
        for (int jj = 0; jj < waveSec.extent(1); ++jj) {
            // factor matrix for calculations
            std::vector<std::vector<std::vector<dcmplx>>> fac;
            fac.assign(rows, std::vector<std::vector<dcmplx>>(cols, std::vector<dcmplx>(depth, dcmplx(0.0, 0.0) ) ) );

            // Define directions for the waves
            double dir1 = dir[ii]; if (materialIndex == 0) dir1 = -1.0;
            double dir2 = dir[jj]; if (sourceIndex == zValsInterfaces.size()) dir2 = -1.0;

            // Distance to interfaces
            auto distZ1 = this->interDist(z1, materialIndex, dir1);
            auto distZ2 = this->interDist(z2, sourceIndex,  -dir2);

            // Add up secondary fields
            for (size_t i = 0; i < rows; i++) {
                for (size_t j = 0; j < cols; j++) {
                    for (size_t k = 0; k < depth; k++) {
                        fac[i][j][k] = I * waveSec(ii, jj) * 
                                       exp(I * (kz1 * distZ1[i][j][k] + kz2 * distZ2[i][j][k]) ) / 
                                       (4 * M_PI * kz2);
                        f[i][j][k]  += fac[i][j][k];
                        
                        // z-derivatives
                        fZ1[i][j][k] = fZ1[i][j][k] + I * kz1 * dir1 * fac[i][j][k];
                        fZ2[i][j][k] = fZ2[i][j][k] - I * kz2 * dir2 * fac[i][j][k];
                        fZZ[i][j][k] = fZZ[i][j][k] + kz1 * dir1 * kz2 * dir2 * fac[i][j][k];
                    }
                }
            }
        }
    }

    return {f, fZ1, fZ2, fZZ};
}

// ----------------------------
// Intra-/Inter-layer integrands
// ----------------------------
/**
 * @details Integrand for Green function evaluation, Chew (6.15, 6.16): 
 * Chew, W.C., Tong, M.S., Hu, B. (2009). 
 * Dyadic Green's Function for Layered Media and Integral Equations. 
 * In: Integral Equation Methods for Electromagnetic and Elastic Waves. 
 * Synthesis Lectures on Computational Electromagnetics. Springer, Cham. 
 * https://doi.org/10.1007/978-3-031-01707-0_6.
 */
IntegrandIntra SommerfeldIntegrator::IntraIntegrand(const std::vector<dcmplx>& rq_te, 
                                                    const std::vector<dcmplx>& rq_tm,
                                                    const Data2D& inData2D, dcmplx kr, dcmplx kz, 
                                                    const std::string& mode, const std::string& op) {
    std::vector<std::vector<double>> r, z;

    // Check if the input data is vectorial or tensorial and convert to 2D arrays
    if (inData2D.rMatrix.empty() && inData2D.z1Matrix.empty())
    {
        // Convert rVec, z1Vec to Nx1 2D arrays.
        size_t N = inData2D.rVec.size();
        if (N != inData2D.z1Vec.size()) {
            // Handle error or throw an exception if mismatch
            throw std::runtime_error("rVec and z1Vec must have the same size");
        }

        r.resize(N, std::vector<double>(1));
        z.resize(N, std::vector<double>(1));

        for (size_t i = 0; i < N; ++i) {
            r[i][0] = inData2D.rVec[i];
            z[i][0] = inData2D.z1Vec[i];
        }
    } else if (inData2D.rVec.empty() && inData2D.z1Vec.empty()) {
        r = inData2D.rMatrix;
        z = inData2D.z1Matrix;
    } else {
        // Handle error or throw an exception if mismatch
        throw std::runtime_error("Data2D structure empty or mismatched");
    }
    size_t rows = r.size();
    size_t cols = r[0].size();

    // Initialize output structure
    IntegrandIntra interIntra;

    // Green function elements for K term in the final system
    interIntra.te.resize(rows, std::vector<dcmplx>(cols));
    interIntra.tm.resize(rows, std::vector<dcmplx>(cols));
    // First derivative along z, subtract quasistatic contribution
    interIntra.tez.resize(rows, std::vector<dcmplx>(cols));
    interIntra.tmz.resize(rows, std::vector<dcmplx>(cols));
    // Second derivative along z
    interIntra.tezz.resize(rows, std::vector<dcmplx>(cols));
    interIntra.tmzz.resize(rows, std::vector<dcmplx>(cols));
    // Surface Green function, Chew (6.28)
    interIntra.tes.resize(rows, std::vector<dcmplx>(cols));
    interIntra.tms.resize(rows, std::vector<dcmplx>(cols));
    // Green function elements for D term in the final system, radial derivative
    interIntra.ter.resize(rows, std::vector<dcmplx>(cols));
    interIntra.tmr.resize(rows, std::vector<dcmplx>(cols));
    // Mixed derivatives
    interIntra.terz.resize(rows, std::vector<dcmplx>(cols));
    interIntra.tmrz.resize(rows, std::vector<dcmplx>(cols));

    // Compute Bessel or Hankel functions
    std::vector<std::vector<dcmplx>> z0, z1;
    z0.resize(rows, std::vector<dcmplx>(cols));                                                                    
    z1.resize(rows, std::vector<dcmplx>(cols));                                                                    

    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            if (mode == "bessel") {
                z0[i][j] = besselJ(0, kr * r[i][j]);
                z1[i][j] = besselJ(1, kr * r[i][j]);
            } else if (mode == "hankel") {
                z0[i][j] = hankelH(0, 1, kr * r[i][j]);
                z1[i][j] = hankelH(1, 1, kr * r[i][j]);
            }
        }
    }

    if (op == "") {
        // Compute Green function elements
        std::vector<std::vector<dcmplx>> f0(rows, std::vector<dcmplx>(cols));
        std::vector<std::vector<dcmplx>> f1(rows, std::vector<dcmplx>(cols));
        std::vector<std::vector<dcmplx>> f0q(rows, std::vector<dcmplx>(cols));
        std::vector<std::vector<dcmplx>> f1q(rows, std::vector<dcmplx>(cols));

        dcmplx q = sqrt(-kr * kr);          // For quasistatic approximation
        if (std::imag(q) < 0)  q = -q;      // Ensure positive imaginary part
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                // Integrand w/o reflection coefficients
                f0[i][j] = I * exp(I * kz * z[i][j]) / (4.0 * M_PI * kz) * z0[i][j];
                f1[i][j] = I * exp(I * kz * z[i][j]) / (4.0 * M_PI * kz) * z1[i][j];
                // Quasistatic approximation
                f0q[i][j] = I * exp(I * q * z[i][j]) / (4.0 * M_PI) * z0[i][j];
                f1q[i][j] = I * exp(I * q * z[i][j]) / (4.0 * M_PI) * z1[i][j] * 
                            double(this->layeredUtils->sgn<double>(std::real(kr)));
            }
        }
        
        // Compute generalized reflection coefficients
        dcmplx r_te, r_tm;
        double dir = 0.0;
        if (materialIndex == 0) {
            std::vector<dcmplx> teRT(this->layeredUtils->totalFresnelCoeffs(kr, "TE", "up"));
            std::vector<dcmplx> tmRT(this->layeredUtils->totalFresnelCoeffs(kr, "TM", "up"));
            r_te = teRT[0]; r_tm = tmRT[0]; dir = -1.0;
        } else {
            std::vector<dcmplx> teRT(this->layeredUtils->totalFresnelCoeffs(kr, "TE", "down"));
            std::vector<dcmplx> tmRT(this->layeredUtils->totalFresnelCoeffs(kr, "TM", "down"));
            r_te = teRT[0]; r_tm = tmRT[0]; dir = 1.0;
        }
        
        // Compute Green function elements
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                // Green function elements for K term in the final system
                interIntra.te[i][j] = r_te * f0[i][j];
                interIntra.tm[i][j] = r_tm * f0[i][j];
                // First derivative along z, subtract quasistatic contribution
                interIntra.tez[i][j] = dir * (I * kz * r_te * f0[i][j] - I * rq_te[0] * f0q[i][j]);
                interIntra.tmz[i][j] = dir * (I * kz * r_tm * f0[i][j] - I * rq_tm[0] * f0q[i][j]);
                // Second derivative along z
                interIntra.tezz[i][j] = -kz * kz * r_te * f0[i][j] + q * rq_te[0] * f0q[i][j];
                interIntra.tmzz[i][j] = -kz * kz * r_tm * f0[i][j] + q * rq_tm[0] * f0q[i][j];
                // Surface Green function, Chew (6.28)
                interIntra.tes[i][j] = kr * kr * r_te * f0[i][j] + q * rq_te[0] * f0q[i][j];
                interIntra.tms[i][j] = kr * kr * r_tm * f0[i][j] + q * rq_tm[0] * f0q[i][j];
                // Green function elements for D term in the final system, radial derivative
                interIntra.ter[i][j] = -kr * r_te * f1[i][j] - I * rq_te[0] * f1q[i][j];
                interIntra.tmr[i][j] = -kr * r_tm * f1[i][j] - I * rq_tm[0] * f1q[i][j];
                // Mixed derivatives
                interIntra.terz[i][j] = dir * (-I * kz * kr * r_te * f1[i][j] + q * rq_te[0] * f1q[i][j]);
                interIntra.tmrz[i][j] = dir * (-I * kz * kr * r_tm * f1[i][j] + q * rq_tm[0] * f1q[i][j]);
            }
        }
    } else {
        // Layer thickness
        double thick = thickness[this->materialIndex - 1];

        // Compute Green function elements
        std::vector<std::vector<dcmplx>> f1(rows, std::vector<dcmplx>(cols));
        std::vector<std::vector<dcmplx>> f2(rows, std::vector<dcmplx>(cols));
        std::vector<std::vector<dcmplx>> f1q(rows, std::vector<dcmplx>(cols));
        std::vector<std::vector<dcmplx>> f2q(rows, std::vector<dcmplx>(cols));

        dcmplx q = sqrt(-kr * kr);          // For quasistatic approximation
        if (std::imag(q) < 0)  q = -q;      // Ensure positive imaginary part
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                // Up- and downgoing waves
                f1[i][j] = I * exp(I * kz * z[i][j]) / (4.0 * M_PI * kz);
                f2[i][j] = I * exp(I * kz * (2.0 *  thick - z[i][j])) / (4.0 * M_PI * kz);
                // Quasistatic approximation
                f1q[i][j] = I * exp(I * q * z[i][j]) / (4.0 * M_PI);
                f2q[i][j] = I * exp(I * q * (2.0 *  thick - z[i][j])) / (4.0 * M_PI);
            }
        }
        double s = double(this->layeredUtils->sgn<double>(std::real(kr)));

        // Secondary waves
        blitz::Array<dcmplx, 2> waveSecTE = 
                    this->layeredUtils->Secondary(kr, this->materialIndex, this->materialIndex, "TE");
        blitz::Array<dcmplx, 2> waveSecTM = 
                    this->layeredUtils->Secondary(kr, this->materialIndex, this->materialIndex, "TM");    
        
        if (op == "minus") {
            // Compute Green function elements
            for (size_t i = 0; i < rows; ++i) {
                for (size_t j = 0; j < cols; ++j) {
                    // Green function elements for K term in the final system
                    interIntra.te[i][j] = (waveSecTE(0,0) * f1[i][j] + waveSecTE(1,1) * f2[i][j]) * z0[i][j];
                    interIntra.tm[i][j] = (waveSecTM(0,0) * f1[i][j] + waveSecTM(1,1) * f2[i][j]) * z0[i][j];
                    // First derivative along z, subtract quasistatic contribution
                    interIntra.tez[i][j] = I * kz * (waveSecTE(0,0) * f1[i][j] - waveSecTE(1,1) * f2[i][j]) * z0[i][j];
                    interIntra.tmz[i][j] = I * kz * (waveSecTM(0,0) * f1[i][j] - waveSecTM(1,1) * f2[i][j]) * z0[i][j];
                    // Second derivative along z
                    interIntra.tezz[i][j] = -kz * kz * (waveSecTE(0,0) * f1[i][j] + waveSecTE(1,1) * f2[i][j] ) * z0[i][j];
                    interIntra.tmzz[i][j] = -kz * kz * (waveSecTM(0,0) * f1[i][j] + waveSecTM(1,1) * f2[i][j] ) * z0[i][j];
                    // Surface Green function, Chew (6.28)
                    interIntra.tes[i][j] = kr * kr * (waveSecTE(0,0) * f1[i][j] + waveSecTE(1,1) * f2[i][j] ) * z0[i][j];
                    interIntra.tms[i][j] = kr * kr * (waveSecTM(0,0) * f1[i][j] + waveSecTM(1,1) * f2[i][j] ) * z0[i][j];
                    // Green function elements for D term in the final system, radial derivative
                    interIntra.ter[i][j] = -kr * (waveSecTE(0,0) * f1[i][j] + waveSecTE(1,1) * f2[i][j] ) * z1[i][j];
                    interIntra.tmr[i][j] = -kr * (waveSecTM(0,0) * f1[i][j] + waveSecTM(1,1) * f2[i][j] ) * z1[i][j];
                    // Mixed derivatives
                    interIntra.terz[i][j] = -I * kz * kr * (waveSecTE(0,0) * f1[i][j] - waveSecTE(1,1) * f2[i][j] ) * z1[i][j];
                    interIntra.tmrz[i][j] = -I * kz * kr * (waveSecTM(0,0) * f1[i][j] - waveSecTM(1,1) * f2[i][j] ) * z1[i][j];
                }
            }
        } else if (op == "plus") {
            // Compute Green function elements
            for (size_t i = 0; i < rows; ++i) {
                for (size_t j = 0; j < cols; ++j) {
                    // Green function elements for K term in the final system
                    interIntra.te[i][j] = (waveSecTE(0,1) * f1[i][j] + waveSecTE(1,0) * f2[i][j]) * z0[i][j];
                    interIntra.tm[i][j] = (waveSecTM(0,1) * f1[i][j] + waveSecTM(1,0) * f2[i][j]) * z0[i][j];
                    // First derivative along z, subtract quasistatic contribution
                    interIntra.tez[i][j] = I * kz * (waveSecTE(0,1) * f1[i][j] - waveSecTE(1,0) * f2[i][j]) * z0[i][j] - 
                                            I * (rq_te[0] * f1q[i][j] - rq_te[1] * f2q[i][j]) * z0[i][j];
                    interIntra.tmz[i][j] = I * kz * (waveSecTM(0,1) * f1[i][j] - waveSecTM(1,0) * f2[i][j]) * z0[i][j] - 
                                            I * (rq_tm[0] * f1q[i][j] - rq_tm[1] * f2q[i][j]) * z0[i][j]; 
                    // Second derivative along z
                    interIntra.tezz[i][j] = -kz * kz * (waveSecTE(0,1) * f1[i][j] + waveSecTE(1,0) * f2[i][j] ) * z0[i][j] + 
                                                q * (rq_te[0] * f1q[i][j] + rq_te[1] * f2q[i][j]) * z0[i][j];
                    interIntra.tmzz[i][j] = -kz * kz * (waveSecTM(0,1) * f1[i][j] + waveSecTM(1,0) * f2[i][j] ) * z0[i][j] + 
                                                q * (rq_tm[0] * f1q[i][j] + rq_tm[1] * f2q[i][j]) * z0[i][j];
                    // Surface Green function, Chew (6.28)
                    interIntra.tes[i][j] = kr * kr * (waveSecTE(0,1) * f1[i][j] + waveSecTE(1,0) * f2[i][j] ) * z0[i][j] + 
                                            q * (rq_te[0] * f1q[i][j] + rq_te[1] * f2q[i][j] ) * z0[i][j];
                    interIntra.tms[i][j] = kr * kr * (waveSecTM(0,1) * f1[i][j] + waveSecTM(1,0) * f2[i][j] ) * z0[i][j] + 
                                            q * (rq_tm[0] * f1q[i][j] + rq_tm[1] * f2q[i][j] ) * z0[i][j];
                    // Green function elements for D term in the final system, radial derivative
                    interIntra.ter[i][j] = -kr * (waveSecTE(0,1) * f1[i][j] + waveSecTE(1,0) * f2[i][j] ) * z1[i][j] - 
                                                I * s * (rq_te[0] * f1q[i][j] + rq_te[1] * f2q[i][j] ) * z1[i][j];
                    interIntra.tmr[i][j] = -kr * (waveSecTM(0,1) * f1[i][j] + waveSecTM(1,0) * f2[i][j] ) * z1[i][j] - 
                                                I * s * (rq_tm[0] * f1q[i][j] + rq_tm[1] * f2q[i][j] ) * z1[i][j];
                    // Mixed derivatives
                    interIntra.terz[i][j] = -I * kz * kr * (waveSecTE(0,1) * f1[i][j] - waveSecTE(1,0) * f2[i][j] ) * z1[i][j] + 
                                                q * s * (rq_te[0] * f1q[i][j] - rq_te[1] * f2q[i][j] ) * z1[i][j];
                    interIntra.tmrz[i][j] = -I * kz * kr * (waveSecTM(0,1) * f1[i][j] - waveSecTM(1,0) * f2[i][j] ) * z1[i][j] + 
                                                q * s * (rq_tm[0] * f1q[i][j] - rq_tm[1] * f2q[i][j] ) * z1[i][j];
                }
            }
        } else {
            // Handle error or throw an exception if mismatch
            throw std::runtime_error("Wrong type for intra-layer Green's function calcularions: use 'plus' or 'minus'.");
        }
    }

    return interIntra;
}

/**
 * @details Integrand for Green function evaluation, Chew (6.15, 6.16): 
 * Chew, W.C., Tong, M.S., Hu, B. (2009). 
 * Dyadic Green's Function for Layered Media and Integral Equations. 
 * In: Integral Equation Methods for Electromagnetic and Elastic Waves. 
 * Synthesis Lectures on Computational Electromagnetics. Springer, Cham. 
 * https://doi.org/10.1007/978-3-031-01707-0_6.
 */
IntegrandInter SommerfeldIntegrator::InterIntegrand(const dcmplx& tq_te, const dcmplx& tq_tm,
                                                    const Data3D& inData3D, dcmplx kr, 
                                                    const std::string& mode) {
    std::vector<std::vector<std::vector<double>>> r, z1, z2;

    // Check if the input data is vectorial or tensorial and convert to 2D arrays
    if (inData3D.rMatrix.empty() && inData3D.z1Matrix.empty() && inData3D.z2Matrix.empty())
    {
        // Convert rVec, z1Vec to Nx1 2D arrays.
        size_t N = inData3D.rVec.size();
        if (N != inData3D.z1Vec.size() || N != inData3D.z2Vec.size()) {
            // Handle error or throw an exception if mismatch
            throw std::runtime_error("rVec and z1Vec must have the same size");
        }

        r.resize( N, std::vector<std::vector<double>>(1, std::vector<double>(1)));
        z1.resize(N, std::vector<std::vector<double>>(1, std::vector<double>(1)));
        z2.resize(N, std::vector<std::vector<double>>(1, std::vector<double>(1)));

        for (size_t i = 0; i < N; ++i) {
            r[i][0][0]  = inData3D.rVec[i];
            z1[i][0][0] = inData3D.z1Vec[i];
            z2[i][0][0] = inData3D.z2Vec[i];
        }
    } else if (inData3D.rVec.empty() && inData3D.z1Vec.empty() && inData3D.z2Vec.empty()) {
        r  = inData3D.rMatrix;
        z1 = inData3D.z1Matrix;
        z2 = inData3D.z2Matrix;
    } else {
        // Handle error or throw an exception if mismatch
        throw std::runtime_error("Data3D structure empty or mismatched");
    }
    size_t rows  = r.size();
    size_t cols  = r[0].size();
    size_t depth = r[0][0].size();

    // Initialize output structure
    IntegrandInter interInter;

    // Green function elements for K term in the final system
    interInter.te.resize(rows, std::vector<std::vector<dcmplx>>(cols, std::vector<dcmplx>(depth)));
    interInter.tm.resize(rows, std::vector<std::vector<dcmplx>>(cols, std::vector<dcmplx>(depth)));
    // First derivative along z, subtract quasistatic contribution
    interInter.tez1.resize(rows, std::vector<std::vector<dcmplx>>(cols, std::vector<dcmplx>(depth)));
    interInter.tez2.resize(rows, std::vector<std::vector<dcmplx>>(cols, std::vector<dcmplx>(depth)));
    interInter.tmz1.resize(rows, std::vector<std::vector<dcmplx>>(cols, std::vector<dcmplx>(depth)));
    interInter.tmz2.resize(rows, std::vector<std::vector<dcmplx>>(cols, std::vector<dcmplx>(depth)));
    // Second derivative along z
    interInter.tezz.resize(rows, std::vector<std::vector<dcmplx>>(cols, std::vector<dcmplx>(depth)));
    interInter.tmzz.resize(rows, std::vector<std::vector<dcmplx>>(cols, std::vector<dcmplx>(depth)));
    // Surface Green function, Chew (6.28)
    interInter.tes.resize(rows, std::vector<std::vector<dcmplx>>(cols, std::vector<dcmplx>(depth)));
    interInter.tms.resize(rows, std::vector<std::vector<dcmplx>>(cols, std::vector<dcmplx>(depth)));
    // Green function elements for D term in the final system, radial derivative
    interInter.ter.resize(rows, std::vector<std::vector<dcmplx>>(cols, std::vector<dcmplx>(depth)));
    interInter.tmr.resize(rows, std::vector<std::vector<dcmplx>>(cols, std::vector<dcmplx>(depth)));
    // Mixed derivatives
    interInter.terz1.resize(rows, std::vector<std::vector<dcmplx>>(cols, std::vector<dcmplx>(depth)));
    interInter.terz2.resize(rows, std::vector<std::vector<dcmplx>>(cols, std::vector<dcmplx>(depth)));
    interInter.tmrz1.resize(rows, std::vector<std::vector<dcmplx>>(cols, std::vector<dcmplx>(depth)));
    interInter.tmrz2.resize(rows, std::vector<std::vector<dcmplx>>(cols, std::vector<dcmplx>(depth)));

    // Compute Bessel or Hankel functions
    std::vector<std::vector<std::vector<dcmplx>>> b0, b1;
    b0.resize(rows, std::vector<std::vector<dcmplx>>(cols, std::vector<dcmplx>(depth)));                                                                    
    b1.resize(rows, std::vector<std::vector<dcmplx>>(cols, std::vector<dcmplx>(depth)));                                                                    

    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            for (size_t k = 0; k < depth; ++k) {
                if (mode == "bessel") {
                    b0[i][j][k] = besselJ(0, kr * r[i][j][k]);
                    b1[i][j][k] = besselJ(1, kr * r[i][j][k]);
                } else if (mode == "hankel") {
                    b0[i][j][k] = hankelH(0, 1, kr * r[i][j][k]);
                    b1[i][j][k] = hankelH(1, 1, kr * r[i][j][k]);
                }
            }
        }
    }

    std::vector<std::vector<std::vector<std::vector<dcmplx>>>> fAll;
    fAll = this->Secondary(kr, r, z1, z2, "TE");
    std::vector<std::vector<std::vector<dcmplx>>> fTE   = fAll[0], fTEZ1 = fAll[1];
    std::vector<std::vector<std::vector<dcmplx>>> fTEZ2 = fAll[2], fTEZZ = fAll[3];
    fAll = this->Secondary(kr, r, z1, z2, "TM");
    std::vector<std::vector<std::vector<dcmplx>>> fTM   = fAll[0], fTMZ1 = fAll[1];
    std::vector<std::vector<std::vector<dcmplx>>> fTMZ2 = fAll[2], fTMZZ = fAll[3];

    // Propagation direction of reflected waves and propagation distance
    double dir = this->layeredUtils->sgn<double>(double(materialIndex - sourceIndex));
    std::vector<std::vector<std::vector<double>>> Z;
    Z.resize(rows, std::vector<std::vector<double>>(cols, std::vector<double>(depth)));
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            for (size_t k = 0; k < depth; ++k) {
                Z[i][j][k] = std::abs(z1[i][j][k] - z2[i][j][k]);
            }
        }
    }

    // Compute Green function elements
    std::vector<std::vector<std::vector<dcmplx>>> f0q, f1q;
    f0q.resize(rows, std::vector<std::vector<dcmplx>>(cols, std::vector<dcmplx>(depth)));
    f1q.resize(rows, std::vector<std::vector<dcmplx>>(cols, std::vector<dcmplx>(depth)));

    dcmplx q = sqrt(-kr * kr);          // For quasistatic approximation
    if (std::imag(q) < 0)  q = -q;      // Ensure positive imaginary part
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            for (size_t k = 0; k < depth; ++k) {
                // Up- and downgoing waves
                f0q[i][j][k] = I * exp(I * q * Z[i][j][k]) / (4.0 * M_PI) * b0[i][j][k];
                f1q[i][j][k] = I * exp(I * q * Z[i][j][k]) / (4.0 * M_PI) * b1[i][j][k] * 
                               ( double(this->layeredUtils->sgn<double>( std::real(kr) ) ) );
            }
        }
    }
    // Compute Green function elements
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            for (size_t k = 0; k < depth; ++k) {
                // Green function elements for K term in the final system
                interInter.te[i][j][k] = fTE[i][j][k] * b0[i][j][k];
                interInter.tm[i][j][k] = fTM[i][j][k] * b0[i][j][k];
                // First derivative along z, subtract quasistatic contribution
                interInter.tez1[i][j][k] = fTEZ1[i][j][k] * b0[i][j][k] - I * dir * tq_te * f0q[i][j][k];
                interInter.tmz1[i][j][k] = fTMZ1[i][j][k] * b0[i][j][k] - I * dir * tq_tm * f0q[i][j][k];
                interInter.tez2[i][j][k] = fTEZ2[i][j][k] * b0[i][j][k] + I * dir * tq_te * f0q[i][j][k];
                interInter.tmz2[i][j][k] = fTMZ2[i][j][k] * b0[i][j][k] + I * dir * tq_tm * f0q[i][j][k];
                // Second derivative along z
                interInter.tezz[i][j][k] = fTEZZ[i][j][k] * b0[i][j][k] - q * tq_te * f0q[i][j][k];
                interInter.tmzz[i][j][k] = fTMZZ[i][j][k] * b0[i][j][k] - q * tq_tm * f0q[i][j][k];
                // Surface Green function, Chew (6.28)
                interInter.tes[i][j][k] = kr * kr * fTE[i][j][k] * b0[i][j][k] + q * tq_te * f0q[i][j][k];
                interInter.tms[i][j][k] = kr * kr * fTM[i][j][k] * b0[i][j][k] + q * tq_tm * f0q[i][j][k];
                // Green function elements for D term in the final system, radial derivative
                interInter.ter[i][j][k] = -kr * fTE[i][j][k] * b1[i][j][k] - I * tq_te * f1q[i][j][k];
                interInter.tmr[i][j][k] = -kr * fTM[i][j][k] * b1[i][j][k] - I * tq_tm * f1q[i][j][k];
                // Mixed derivatives
                interInter.terz1[i][j][k] = -kr * fTEZ1[i][j][k] * b1[i][j][k] + dir * q * tq_te * f1q[i][j][k];
                interInter.tmrz1[i][j][k] = -kr * fTMZ1[i][j][k] * b1[i][j][k] + dir * q * tq_tm * f1q[i][j][k];
                interInter.terz2[i][j][k] = -kr * fTEZ2[i][j][k] * b1[i][j][k] - dir * q * tq_te * f1q[i][j][k];
                interInter.tmrz2[i][j][k] = -kr * fTMZ2[i][j][k] * b1[i][j][k] - dir * q * tq_tm * f1q[i][j][k];
            }
        }
    }
    return interInter;
}

// ----------------------------
// One-abscissa drivers w/ subtraction
// ----------------------------
// Integration function for internal intra-layer Sommerfeld integrals.
IntegrandIntra SommerfeldIntegrator::integrationFIntra(const std::vector<dcmplx>& rq_te, 
                                                       const std::vector<dcmplx>& rq_tm, 
                                                       const Data2D& inData, double x, double k_max,
                                                       const IntegrandIntra& intZeroKrho, 
                                                       const std::string& path, const std::string& op) {
    IntegrandIntra result;
    dcmplx k_layer = this->k_L[materialIndex];
    double a = this->cutoff * 
               std::real( k_layer / ( sqrt( this->epsilon[materialIndex] * this->mu[materialIndex] ) ) );
    
    // Integration along semi-ellipse
    if (path == "semi") {
        // krho parametric form
        dcmplx kr = k_max * (1.0 - cos(x) - I * this->semiImag * sin(x));

        // krho parametric form derivative w.r.t. x
        dcmplx dkr_dx = k_max * (sin(x) - I * this->semiImag * cos(x));

        // kz of the layer
        dcmplx kz = sqrt(k_layer * k_layer - kr * kr);
        if (std::imag(kz) < 0) kz = -kz;

        // Compute the integrand
        IntegrandIntra intSemi = IntraIntegrand(rq_te, rq_tm, inData, kr, kz, "bessel", op);
        dcmplx factorSemi = dkr_dx / kr;
        intSemi = factorSemi * intSemi;

        // Singularity subtraction
        factorSemi = dkr_dx / kr * a * a / (a * a + kr * kr);
        result = intSemi - factorSemi * intZeroKrho;
    }
    // Integration to real infinity
    else if (path == "real") {
        // krho parametric form
        dcmplx kr = 2.0 * k_max / x;

        // krho parametric form derivative w.r.t. x
        dcmplx dkr_dx = -2.0 * k_max / (x * x);

        // kz of the layer
        dcmplx kz = sqrt(k_layer * k_layer - kr * kr);
        if (std::imag(kz) < 0) kz = -kz;

        // Select the corresponding data for the integration path
        IntegrandIntra intZeroKrhoReal = intZeroKrho.flatten(inData.indVecComplex);

        // Compute the integrand
        IntegrandIntra intReal = IntraIntegrand(rq_te, rq_tm, inData, kr, kz, "bessel", op);
        dcmplx factorReal = dkr_dx / kr;
        intReal = factorReal * intReal;
        factorReal = dkr_dx / kr * (a * a) / (a * a + kr * kr);

        // Singularity subtraction
        result = intReal - factorReal * intZeroKrhoReal;
    }
    
    // Integration to imaginary infinity
    else if (path == "imag") {
        // krho parametric form
        dcmplx kr1 = 2.0 * k_max * ( 1.0 - I + I / x);
        dcmplx kr2 = 2.0 * k_max * (-1.0 - I + I / x);

        // krho parametric form derivative w.r.t. x
        dcmplx dkr_dx = -2.0 * I * k_max / (x * x);
        
        // kz of the layer
        dcmplx kz1 = sqrt(k_layer * k_layer - kr1 * kr1); if (std::imag(kz1) < 0) kz1 = -kz1;
        dcmplx kz2 = sqrt(k_layer * k_layer - kr2 * kr2); if (std::imag(kz2) < 0) kz2 = -kz2;
        
        // Select the corresponding data for the integration path
        IntegrandIntra intZeroKrhoImag = intZeroKrho.flatten(inData.indVecComplex);

        // Compute the integrands
        IntegrandIntra intImag1 = IntraIntegrand(rq_te, rq_tm, inData, kr1, kz1, "hankel", op);
        IntegrandIntra intImag2 = IntraIntegrand(rq_te, rq_tm, inData, kr2, kz2, "hankel", op);
        dcmplx factorImag1 = 0.5 * dkr_dx / kr1;
        dcmplx factorImag2 = 0.5 * dkr_dx / kr2;

        IntegrandIntra intImag = factorImag1 * intImag1 - factorImag2 * intImag2;
        factorImag1 = dkr_dx / kr1 * (a * a) / (a * a + kr1 * kr1);
        factorImag2 = dkr_dx / kr2 * (a * a) / (a * a + kr2 * kr2);
        
        // Singularity subtraction
        result = intImag - (factorImag1 - factorImag2) * intZeroKrhoImag;
    }
    // Invalid path type
    else {
        throw std::invalid_argument("Invalid path type. Use 'semi', 'real', or 'imag'.");
    }
    
    return result;
}

// Integration function for inter-layer Sommerfeld integrals.
IntegrandInter SommerfeldIntegrator::integrationFInter(const dcmplx& tq_te, const dcmplx& tq_tm, 
                                                       const Data3D& inData, double x, double k_max,
                                                       const IntegrandInter& intZeroKrho, 
                                                       const std::string& path) {
    IntegrandInter result;
    dcmplx k_layer = this->k_L[materialIndex];
    double a = this->cutoff * 
               std::real( k_layer / ( sqrt( this->epsilon[materialIndex] * this->mu[materialIndex] ) ) );
    
    // Integration along semi-ellipse
    if (path == "semi") {
        // krho parametric form
        dcmplx kr = k_max * (1.0 - cos(x) - I * this->semiImag * sin(x));

        // krho parametric form derivative w.r.t. x
        dcmplx dkr_dx = k_max * (sin(x) - I * this->semiImag * cos(x));

        // Compute the integrand
        IntegrandInter intSemi = InterIntegrand(tq_te, tq_tm, inData, kr, "bessel");
        dcmplx factorSemi = dkr_dx / kr;
        intSemi = factorSemi * intSemi;

        // Singularity subtraction
        factorSemi = dkr_dx / kr * a * a / (a * a + kr * kr);
        result = intSemi - factorSemi * intZeroKrho;
    }
    // Integration to real infinity
    else if (path == "real") {
        // krho parametric form
        dcmplx kr = 2.0 * k_max / x;

        // krho parametric form derivative w.r.t. x
        dcmplx dkr_dx = -2.0 * k_max / (x * x);

        // Select the corresponding data for the integration path
        IntegrandInter intZeroKrhoReal = intZeroKrho.flatten(inData.indVecTriplet);

        // Compute the integrand
        IntegrandInter intReal = InterIntegrand(tq_te, tq_tm, inData, kr, "bessel");
        dcmplx factorReal = dkr_dx / kr;
        intReal = factorReal * intReal;
        factorReal = dkr_dx / kr * (a * a) / (a * a + kr * kr);

        // Singularity subtraction
        result = intReal - factorReal * intZeroKrhoReal;
    }
    
    // Integration to imaginary infinity
    else if (path == "imag") {
        // krho parametric form
        dcmplx kr1 = 2.0 * k_max * ( 1.0 - I + I / x);
        dcmplx kr2 = 2.0 * k_max * (-1.0 - I + I / x);

        // krho parametric form derivative w.r.t. x
        dcmplx dkr_dx = -2.0 * I * k_max / (x * x);
        
        // Select the corresponding data for the integration path
        IntegrandInter intZeroKrhoImag = intZeroKrho.flatten(inData.indVecTriplet);

        // Compute the integrands
        IntegrandInter intImag1 = InterIntegrand(tq_te, tq_tm, inData, kr1, "hankel");
        IntegrandInter intImag2 = InterIntegrand(tq_te, tq_tm, inData, kr2, "hankel");
        dcmplx factorImag1 = 0.5 * dkr_dx / kr1;
        dcmplx factorImag2 = 0.5 * dkr_dx / kr2;

        IntegrandInter intImag = factorImag1 * intImag1 - factorImag2 * intImag2;
        factorImag1 = dkr_dx / kr1 * (a * a) / (a * a + kr1 * kr1);
        factorImag2 = dkr_dx / kr2 * (a * a) / (a * a + kr2 * kr2);
        
        // Singularity subtraction
        result = intImag - (factorImag1 - factorImag2) * intZeroKrhoImag;
    }
    // Invalid path type
    else {
        throw std::invalid_argument("Invalid path type. Use 'semi', 'real', or 'imag'.");
    }
    
    return result;
}

/**
 * @brief Extract a single column from a 2D matrix.
 *
 * @param matrix Input 2D matrix (matrix[row][col]).
 * @param colIndex Index of the column to extract.
 * @return A vector containing all elements from the specified column.
 * @throws std::out_of_range if @p colIndex is greater than the number of columns.
 */
std::vector<double> getColumn(const std::vector<std::vector<double>>& matrix, size_t colIndex) {
    std::vector<double> column(matrix.size()); // Preallocate space
    std::transform(matrix.begin(), matrix.end(), column.begin(),
                   [colIndex](const std::vector<double>& row) { return row[colIndex]; });
    return column;
}

/**
 * @brief Extract one row (across all columns) at a specific depth in a 3D matrix.
 *
 * @param matrix Input 3D matrix (matrix[row][col][depth]).
 * @param rowIndex Row index to extract.
 * @param depthIndex Depth (third-dimension) index to extract.
 * @return A vector of values at the specified (row, :, depth).
 * @throws std::out_of_range if indices are out of bounds.
 */
std::vector<double> getColumnAtRowDepth(
    const std::vector<std::vector<std::vector<double>>>& matrix, 
    size_t rowIndex, size_t depthIndex) 
{
    if (rowIndex >= matrix.size()) {
        throw std::out_of_range("Row index out of range.");
    }

    std::vector<double> rowAtDepth;
    for (const auto& col : matrix[rowIndex]) {
        if (depthIndex >= col.size()) {
            throw std::out_of_range("Depth index out of range.");
        }
        rowAtDepth.push_back(col[depthIndex]); // Extract the element at specified depth
    }
    return rowAtDepth;
}

/**
 * @brief Extract all depth elements at a specific row and column in a 3D matrix.
 *
 * @param matrix Input 3D matrix (matrix[row][col][depth]).
 * @param rowIndex Row index to extract.
 * @param colIndex Column index to extract.
 * @return A vector containing all depth values at (rowIndex, colIndex, :).
 * @throws std::out_of_range if indices are out of bounds.
 */
std::vector<double> getDepthAtRowColumn(
    const std::vector<std::vector<std::vector<double>>>& matrix, 
    size_t rowIndex, size_t colIndex) 
{
    if (rowIndex >= matrix.size() || colIndex >= matrix[rowIndex].size()) {
        throw std::out_of_range("Row or column index out of range.");
    }

    return matrix[rowIndex][colIndex]; // Extract the entire depth slice
}

/**
 * @brief Extract one column (across all rows) at a specific depth in a 3D matrix.
 *
 * @param matrix Input 3D matrix (matrix[row][col][depth]).
 * @param colIndex Column index to extract.
 * @param depthIndex Depth index to extract.
 * @return A vector containing the values at (:, colIndex, depthIndex).
 * @throws std::out_of_range if indices are out of bounds.
 */
std::vector<double> getRowAtColumnDepth(
    const std::vector<std::vector<std::vector<double>>>& matrix, 
    size_t colIndex, size_t depthIndex) 
{
    std::vector<double> columnAtDepth;
    for (const auto& row : matrix) {
        if (colIndex >= row.size() || depthIndex >= row[colIndex].size()) {
            throw std::out_of_range("Column or depth index out of range.");
        }
        columnAtDepth.push_back(row[colIndex][depthIndex]); // Extract element
    }
    return columnAtDepth;
}

// Sommerfeld integral evaluation methods for (all) the intra-layer integrals.
CoeffsIntra SommerfeldIntegrator::Evaluate(const std::vector<dcmplx>& rq_te, 
                                           const std::vector<dcmplx>& rq_tm,
                                           const std::vector<std::vector<double>>& r, 
                                           const std::vector<std::vector<double>>& z, 
                                           bool singular, const std::string& op) {
    // Validate singular parameter
    if (singular == false) {
        throw std::invalid_argument("Invalid value for singular parameter. Must be true.");
    } else {
        // Get the layers wavenumbers
        std::vector<dcmplx> k = k_L;
        
        // Material index layer wavenumber
        dcmplx k_layer = k[materialIndex];

        // Maximum layer real part of wavenumbers
        dcmplx temp = *std::max_element(k.begin(), k.end(), [](const dcmplx& a, const dcmplx& b) {
            return std::real(a) < std::real(b);
        });
        double k_max = this->semiReal * std::real(temp);

        // Select data for integration paths of semiellipse, and real or imaginary infinity
        auto dataSemi2D = selectData2D(r, z, {}, "semi");
        auto dataReal2D = selectData2D(r, z, {}, "real");
        auto dataImag2D = selectData2D(r, z, {}, "imag");

        // Integrand for krho = 0
        IntegrandIntra intZeroKrho = IntraIntegrand(rq_te, rq_tm, dataSemi2D, 
                                                    dcmplx(0.0, 0.0), k_layer, "bessel", op);

        // alpha parameter for singularity subtraction Chew (Ch.6, A-7)
        double alpha = this->cutoff * std::real( k_layer / ( sqrt( this->epsilon[materialIndex] * 
                                                                   this->mu[materialIndex] ) ) );
        
        // Set the global pointers/variables needed by the function:
        g_integ = this;  g_rq_te = rq_te; g_rq_tm = rq_tm;
        g_status = "intra"; g_k_max = k_max; g_op = op; 
        g_intZeroKrho = intZeroKrho;
        
        // Numerical integration structure to define the size of the numerical integration results and store them
        IntegrandIntra numIntSemi = integrationFIntra(rq_te, rq_tm, dataSemi2D, 0.5, k_max, intZeroKrho, "semi", op);

        // Semi-ellipse integration path
        g_memo.clear();
        g_Data2D = dataSemi2D;
        g_path = "semi";

        // Loop over the members of the IntegrandIntra structure
        for (IntraMember mem : { IntraMember::TE, IntraMember::TM, IntraMember::TEZ, IntraMember::TMZ,
                                 IntraMember::TEZZ, IntraMember::TMZZ, IntraMember::TES, IntraMember::TMS,
                                 IntraMember::TER, IntraMember::TMR, IntraMember::TERZ, IntraMember::TMRZ }) 
        {
            g_currentMember = mem;
            for (size_t i = 0; i < numIntSemi.te.size(); i++) {
                g_currentI = i;
                for (size_t j = 0; j < numIntSemi.te[0].size(); j++) {
                    g_currentJ = j;
                    double realPart = std::real(gaussKronrod(pickRealMember, 0.0, M_PI));
                    double imagPart = std::real(gaussKronrod(pickImagMember, 0.0, M_PI));
                    dcmplx totalVal(realPart, imagPart);
                    switch (g_currentMember)
                    {
                        case IntraMember::TE:
                            numIntSemi.te[i][j]   = totalVal + intZeroKrho.te[i][j]   * std::log(alpha); break;
                        case IntraMember::TM:
                            numIntSemi.tm[i][j]   = totalVal + intZeroKrho.tm[i][j]   * std::log(alpha); break;
                        case IntraMember::TEZ:
                            numIntSemi.tez[i][j]  = totalVal + intZeroKrho.tez[i][j]  * std::log(alpha); break;
                        case IntraMember::TMZ:
                            numIntSemi.tmz[i][j]  = totalVal + intZeroKrho.tmz[i][j]  * std::log(alpha); break;
                        case IntraMember::TEZZ: 
                            numIntSemi.tezz[i][j] = totalVal + intZeroKrho.tezz[i][j] * std::log(alpha); break;
                        case IntraMember::TMZZ:
                            numIntSemi.tmzz[i][j] = totalVal + intZeroKrho.tmzz[i][j] * std::log(alpha); break;
                        case IntraMember::TES: 
                            numIntSemi.tes[i][j]  = totalVal + intZeroKrho.tes[i][j]  * std::log(alpha); break;
                        case IntraMember::TMS: 
                            numIntSemi.tms[i][j]  = totalVal + intZeroKrho.tms[i][j]  * std::log(alpha); break;
                        case IntraMember::TER:  
                            numIntSemi.ter[i][j]  = totalVal + intZeroKrho.ter[i][j]  * std::log(alpha); break;
                        case IntraMember::TMR:  
                            numIntSemi.tmr[i][j]  = totalVal + intZeroKrho.tmr[i][j]  * std::log(alpha); break;
                        case IntraMember::TERZ: 
                            numIntSemi.terz[i][j] = totalVal + intZeroKrho.terz[i][j] * std::log(alpha); break;
                        case IntraMember::TMRZ: 
                            numIntSemi.tmrz[i][j] = totalVal + intZeroKrho.tmrz[i][j] * std::log(alpha); break;
                    }            
                }
            }
        }

        // Real axis integration path
        g_memo.clear();
        g_Data2D = dataReal2D;
        g_path = "real";
        if (!g_Data2D.indVecComplex.empty()) {
            IntegrandIntra intZeroKrhoReal = numIntSemi.flatten(g_Data2D.indVecComplex);
            // Loop over the members of the IntegrandIntra structure
            for (IntraMember mem : { IntraMember::TE, IntraMember::TM, IntraMember::TEZ, IntraMember::TMZ,
                                    IntraMember::TEZZ, IntraMember::TMZZ, IntraMember::TES, IntraMember::TMS,
                                    IntraMember::TER, IntraMember::TMR, IntraMember::TERZ, IntraMember::TMRZ }) {
                g_currentMember = mem;
                for (size_t i = 0; i < g_Data2D.indVecComplex.size(); i++) {
                    g_currentI = i;
                    size_t j = 0;
                    g_currentJ = j;
                    double realPart = std::real(gaussKronrod(pickRealMember, 1.0, 0.0));
                    double imagPart = std::real(gaussKronrod(pickImagMember, 1.0, 0.0));
                    dcmplx totalVal(realPart, imagPart);
                    int row = std::real(g_Data2D.indVecComplex[i]);
                    int col = std::imag(g_Data2D.indVecComplex[i]);
                    switch (g_currentMember)
                    {
                        case IntraMember::TE:
                            numIntSemi.te[row][col]   = totalVal + intZeroKrhoReal.te[i][j]  ; break;
                        case IntraMember::TM:
                            numIntSemi.tm[row][col]   = totalVal + intZeroKrhoReal.tm[i][j]  ; break;
                        case IntraMember::TEZ:
                            numIntSemi.tez[row][col]  = totalVal + intZeroKrhoReal.tez[i][j] ; break;
                        case IntraMember::TMZ:
                            numIntSemi.tmz[row][col]  = totalVal + intZeroKrhoReal.tmz[i][j] ; break;
                        case IntraMember::TEZZ: 
                            numIntSemi.tezz[row][col] = totalVal + intZeroKrhoReal.tezz[i][j]; break;
                        case IntraMember::TMZZ:
                            numIntSemi.tmzz[row][col] = totalVal + intZeroKrhoReal.tmzz[i][j]; break;
                        case IntraMember::TES: 
                            numIntSemi.tes[row][col]  = totalVal + intZeroKrhoReal.tes[i][j] ; break;
                        case IntraMember::TMS: 
                            numIntSemi.tms[row][col]  = totalVal + intZeroKrhoReal.tms[i][j] ; break;
                        case IntraMember::TER:  
                            numIntSemi.ter[row][col]  = totalVal + intZeroKrhoReal.ter[i][j] ; break;
                        case IntraMember::TMR:  
                            numIntSemi.tmr[row][col]  = totalVal + intZeroKrhoReal.tmr[i][j] ; break;
                        case IntraMember::TERZ: 
                            numIntSemi.terz[row][col] = totalVal + intZeroKrhoReal.terz[i][j]; break;
                        case IntraMember::TMRZ: 
                            numIntSemi.tmrz[row][col] = totalVal + intZeroKrhoReal.tmrz[i][j]; break;
                    }
                }
            }
        }

        // Imaginary axis integration path
        g_memo.clear();
        g_Data2D = dataImag2D;
        g_path = "imag";
        if (!g_Data2D.indVecComplex.empty()) {
            IntegrandIntra intZeroKrhoImag = numIntSemi.flatten(g_Data2D.indVecComplex);
            // Loop over the members of the IntegrandIntra structure
            for (IntraMember mem : { IntraMember::TE, IntraMember::TM, IntraMember::TEZ, IntraMember::TMZ,
                                    IntraMember::TEZZ, IntraMember::TMZZ, IntraMember::TES, IntraMember::TMS,
                                    IntraMember::TER, IntraMember::TMR, IntraMember::TERZ, IntraMember::TMRZ }) {
                g_currentMember = mem;
                for (size_t i = 0; i < g_Data2D.indVecComplex.size(); i++) {
                    g_currentI = i;
                    size_t j = 0;
                    g_currentJ = j;
                    double realPart = std::real(gaussKronrod(pickRealMember, 1.0, 0.0));
                    double imagPart = std::real(gaussKronrod(pickImagMember, 1.0, 0.0));
                    dcmplx totalVal(realPart, imagPart);
                    int row = std::real(g_Data2D.indVecComplex[i]);
                    int col = std::imag(g_Data2D.indVecComplex[i]);
                    switch (g_currentMember)
                    {
                        case IntraMember::TE:
                            numIntSemi.te[row][col]   = totalVal + intZeroKrhoImag.te[i][j]  ; break;
                        case IntraMember::TM:
                            numIntSemi.tm[row][col]   = totalVal + intZeroKrhoImag.tm[i][j]  ; break;
                        case IntraMember::TEZ:
                            numIntSemi.tez[row][col]  = totalVal + intZeroKrhoImag.tez[i][j] ; break;
                        case IntraMember::TMZ:
                            numIntSemi.tmz[row][col]  = totalVal + intZeroKrhoImag.tmz[i][j] ; break;
                        case IntraMember::TEZZ: 
                            numIntSemi.tezz[row][col] = totalVal + intZeroKrhoImag.tezz[i][j]; break;
                        case IntraMember::TMZZ:
                            numIntSemi.tmzz[row][col] = totalVal + intZeroKrhoImag.tmzz[i][j]; break;
                        case IntraMember::TES: 
                            numIntSemi.tes[row][col]  = totalVal + intZeroKrhoImag.tes[i][j] ; break;
                        case IntraMember::TMS: 
                            numIntSemi.tms[row][col]  = totalVal + intZeroKrhoImag.tms[i][j] ; break;
                        case IntraMember::TER:  
                            numIntSemi.ter[row][col]  = totalVal + intZeroKrhoImag.ter[i][j] ; break;
                        case IntraMember::TMR:  
                            numIntSemi.tmr[row][col]  = totalVal + intZeroKrhoImag.tmr[i][j] ; break;
                        case IntraMember::TERZ: 
                            numIntSemi.terz[row][col] = totalVal + intZeroKrhoImag.terz[i][j]; break;
                        case IntraMember::TMRZ: 
                            numIntSemi.tmrz[row][col] = totalVal + intZeroKrhoImag.tmrz[i][j]; break;
                    }
                }
            }
        }
#ifdef DEBUG
        if (materialIndex == 1) {
            writeMatrixToFile(numIntSemi.te,   "te.txt");   writeMatrixToFile(numIntSemi.tm,   "tm.txt");
            writeMatrixToFile(numIntSemi.tez,  "tez.txt");  writeMatrixToFile(numIntSemi.tmz,  "tmz.txt");
            writeMatrixToFile(numIntSemi.tezz, "tezz.txt"); writeMatrixToFile(numIntSemi.tmzz, "tmzz.txt");
            writeMatrixToFile(numIntSemi.tes,  "tes.txt");  writeMatrixToFile(numIntSemi.tms,  "tms.txt");
            writeMatrixToFile(numIntSemi.ter,  "ter.txt");  writeMatrixToFile(numIntSemi.tmr,  "tmr.txt");
            writeMatrixToFile(numIntSemi.terz, "terz.txt"); writeMatrixToFile(numIntSemi.tmrz, "tmrz.txt");
        }
#endif
        std::vector<double> x = getColumn(r, 0);
        std::vector<double> y = z[0];

        CoeffsIntra coeffs;
        coeffs.te   = interpolateIntra(x, y, numIntSemi.te);    coeffs.tm   = interpolateIntra(x, y, numIntSemi.tm);
        coeffs.tez  = interpolateIntra(x, y, numIntSemi.tez);   coeffs.tmz  = interpolateIntra(x, y, numIntSemi.tmz);
        coeffs.tezz = interpolateIntra(x, y, numIntSemi.tezz);  coeffs.tmzz = interpolateIntra(x, y, numIntSemi.tmzz);
        coeffs.tes  = interpolateIntra(x, y, numIntSemi.tes);   coeffs.tms  = interpolateIntra(x, y, numIntSemi.tms);
        coeffs.ter  = interpolateIntra(x, y, numIntSemi.ter);   coeffs.tmr  = interpolateIntra(x, y, numIntSemi.tmr);
        coeffs.terz = interpolateIntra(x, y, numIntSemi.terz);  coeffs.tmrz = interpolateIntra(x, y, numIntSemi.tmrz);

        return coeffs;
    }
}

// Sommerfeld integral evaluation methods for the inter-layer integrals.
CoeffsInter SommerfeldIntegrator::Evaluate(const dcmplx& tq_te, const dcmplx& tq_tm,
                                           const std::vector<std::vector<std::vector<double>>>& r, 
                                           const std::vector<std::vector<std::vector<double>>>& z1, 
                                           const std::vector<std::vector<std::vector<double>>>& z2, 
                                           bool singular, const std::string& op) {
    // Validate singular parameter
    if (singular == false) {
        throw std::invalid_argument("Invalid value for singular parameter. Must be true.");
    } else {
        // Get the layers wavenumbers
        std::vector<dcmplx> k = k_L;
        
        // Material index layer wavenumber
        dcmplx k_layer = k[materialIndex];

        // Maximum layer real part of wavenumbers
        dcmplx temp = *std::max_element(k.begin(), k.end(), [](const dcmplx& a, const dcmplx& b) {
            return std::real(a) < std::real(b);
        });
        double k_max = this->semiReal * std::real(temp);

        // Select data for integration paths of semiellipse, and real or imaginary infinity
        auto dataSemi3D = selectData3D(r, z1, z2, "semi");
        auto dataReal3D = selectData3D(r, z1, z2, "real");
        auto dataImag3D = selectData3D(r, z1, z2, "imag");

        // Integrand for krho = 0
        IntegrandInter intZeroKrho = InterIntegrand(tq_te, tq_tm, dataSemi3D, dcmplx(0.0, 0.0), "bessel");

        // alpha parameter for singularity subtraction Chew (Ch.6, A-7)
        double alpha = this->cutoff * std::real( k_layer / ( sqrt( this->epsilon[materialIndex] * 
                                                                   this->mu[materialIndex] ) ) );

        // Set the global pointers/variables needed by the function:
        g_integ = this;     g_tq_te = tq_te; g_tq_tm = tq_tm;
        g_status = "inter"; g_k_max = k_max; g_op = op;
        g_intZeroKrhoInter = intZeroKrho;
       
        // Numerical integration structure to define the size of the numerical integration results and store them
        IntegrandInter numIntSemi = integrationFInter(tq_te, tq_tm, dataSemi3D, 0.5, k_max, intZeroKrho, "semi");

        // Semi-ellipse integration path
        g_memoInter.clear();
        g_Data3D = dataSemi3D;
        g_path = "semi";

        // Before integration loop (allocate expected size)
        if (g_memoInter.empty()) g_memoInter.reserve(numIntSemi.te.size() * 
                                                     numIntSemi.te[0].size() * 
                                                     numIntSemi.te[0][0].size());
        
        // Loop over the members of the IntegrandInter structure
        for (InterMember mem : { InterMember::TE,    InterMember::TM,    InterMember::TEZ1,  InterMember::TMZ1,
                                 InterMember::TEZ2,  InterMember::TMZ2,  InterMember::TEZZ,  InterMember::TMZZ, 
                                 InterMember::TES,   InterMember::TMS,   InterMember::TER,   InterMember::TMR, 
                                 InterMember::TERZ1, InterMember::TMRZ1, InterMember::TERZ2, InterMember::TMRZ2 }) 
        {
            g_currentMemberInter = mem;
            for (size_t i = 0; i < numIntSemi.te.size(); i++) {
                g_currentI = i;
                for (size_t j = 0; j < numIntSemi.te[0].size(); j++) {
                    g_currentJ = j;
                    for (size_t k = 0; k < numIntSemi.te[0][0].size(); k++) {
                        g_currentK = k;
                        double realPart = std::real(gaussKronrod(pickRealMember, 0.0, M_PI, 5, 100, 1.0e-8, 1.0e-8));
                        double imagPart = std::real(gaussKronrod(pickImagMember, 0.0, M_PI, 5, 100, 1.0e-8, 1.0e-8));
                        dcmplx totalVal(realPart, imagPart);
                        switch (g_currentMemberInter)
                        {
                            case InterMember::TE:
                                numIntSemi.te[i][j][k]    = totalVal + intZeroKrho.te[i][j][k]    * std::log(alpha); break;
                            case InterMember::TM:
                                numIntSemi.tm[i][j][k]    = totalVal + intZeroKrho.tm[i][j][k]    * std::log(alpha); break;
                            case InterMember::TEZ1:
                                numIntSemi.tez1[i][j][k]  = totalVal + intZeroKrho.tez1[i][j][k]  * std::log(alpha); break;
                            case InterMember::TMZ1:
                                numIntSemi.tmz1[i][j][k]  = totalVal + intZeroKrho.tmz1[i][j][k]  * std::log(alpha); break;
                            case InterMember::TEZ2:
                                numIntSemi.tez2[i][j][k]  = totalVal + intZeroKrho.tez2[i][j][k]  * std::log(alpha); break;
                            case InterMember::TMZ2:
                                numIntSemi.tmz2[i][j][k]  = totalVal + intZeroKrho.tmz2[i][j][k]  * std::log(alpha); break;
                            case InterMember::TEZZ: 
                                numIntSemi.tezz[i][j][k]  = totalVal + intZeroKrho.tezz[i][j][k]  * std::log(alpha); break;
                            case InterMember::TMZZ:
                                numIntSemi.tmzz[i][j][k]  = totalVal + intZeroKrho.tmzz[i][j][k]  * std::log(alpha); break;
                            case InterMember::TES: 
                                numIntSemi.tes[i][j][k]   = totalVal + intZeroKrho.tes[i][j][k]   * std::log(alpha); break;
                            case InterMember::TMS: 
                                numIntSemi.tms[i][j][k]   = totalVal + intZeroKrho.tms[i][j][k]   * std::log(alpha); break;
                            case InterMember::TER:  
                                numIntSemi.ter[i][j][k]   = totalVal + intZeroKrho.ter[i][j][k]   * std::log(alpha); break;
                            case InterMember::TMR:  
                                numIntSemi.tmr[i][j][k]   = totalVal + intZeroKrho.tmr[i][j][k]   * std::log(alpha); break;
                            case InterMember::TERZ1: 
                                numIntSemi.terz1[i][j][k] = totalVal + intZeroKrho.terz1[i][j][k] * std::log(alpha); break;
                            case InterMember::TMRZ1: 
                                numIntSemi.tmrz1[i][j][k] = totalVal + intZeroKrho.tmrz1[i][j][k] * std::log(alpha); break;
                            case InterMember::TERZ2: 
                                numIntSemi.terz2[i][j][k] = totalVal + intZeroKrho.terz2[i][j][k] * std::log(alpha); break;
                            case InterMember::TMRZ2: 
                                numIntSemi.tmrz2[i][j][k] = totalVal + intZeroKrho.tmrz2[i][j][k] * std::log(alpha); break;
                        }            
                    }
                }
            }
        }
        
        // Real axis integration path
        g_memoInter.clear();
        g_Data3D = dataReal3D;
        g_path = "real";
        if (g_memoInter.empty()) g_memoInter.reserve(g_Data3D.indVecTriplet.size()); // Before integration loop (expected size)
        if (!g_Data3D.indVecTriplet.empty()) {
            IntegrandInter intZeroKrhoReal = numIntSemi.flatten(g_Data3D.indVecTriplet);
            // Loop over the members of the IntegrandInter structure
            for (InterMember mem : { InterMember::TE,    InterMember::TM,    InterMember::TEZ1,  InterMember::TMZ1,
                                     InterMember::TEZ2,  InterMember::TMZ2,  InterMember::TEZZ,  InterMember::TMZZ, 
                                     InterMember::TES,   InterMember::TMS,   InterMember::TER,   InterMember::TMR, 
                                     InterMember::TERZ1, InterMember::TMRZ1, InterMember::TERZ2, InterMember::TMRZ2 }) 
            {
                g_currentMemberInter = mem;
                for (size_t i = 0; i < g_Data3D.indVecTriplet.size(); i++) {
                    g_currentI = i; size_t j = 0; g_currentJ = j; size_t k = 0; g_currentK = k;
                    double realPart = std::real(gaussKronrod(pickRealMember, 1.0, 0.0, 5, 100, 1.0e-8, 1.0e-8));
                    double imagPart = std::real(gaussKronrod(pickImagMember, 1.0, 0.0, 5, 100, 1.0e-8, 1.0e-8));
                    dcmplx totalVal(realPart, imagPart);
                    int row = std::get<0>(g_Data3D.indVecTriplet[i]);
                    int col = std::get<1>(g_Data3D.indVecTriplet[i]);
                    int dep = std::get<2>(g_Data3D.indVecTriplet[i]);
                    switch (g_currentMemberInter)
                    {
                        case InterMember::TE:
                            numIntSemi.te[row][col][dep]    = totalVal + intZeroKrhoReal.te[i][j][k]  ;  break;
                        case InterMember::TM:
                            numIntSemi.tm[row][col][dep]    = totalVal + intZeroKrhoReal.tm[i][j][k]  ;  break;
                        case InterMember::TEZ1:
                            numIntSemi.tez1[row][col][dep]  = totalVal + intZeroKrhoReal.tez1[i][j][k] ; break;
                        case InterMember::TMZ1:
                            numIntSemi.tmz1[row][col][dep]  = totalVal + intZeroKrhoReal.tmz1[i][j][k] ; break;
                        case InterMember::TEZ2:
                            numIntSemi.tez2[row][col][dep]  = totalVal + intZeroKrhoReal.tez2[i][j][k] ; break;
                        case InterMember::TMZ2:
                            numIntSemi.tmz2[row][col][dep]  = totalVal + intZeroKrhoReal.tmz2[i][j][k] ; break;
                        case InterMember::TEZZ: 
                            numIntSemi.tezz[row][col][dep]  = totalVal + intZeroKrhoReal.tezz[i][j][k];  break;
                        case InterMember::TMZZ:
                            numIntSemi.tmzz[row][col][dep]  = totalVal + intZeroKrhoReal.tmzz[i][j][k];  break;
                        case InterMember::TES: 
                            numIntSemi.tes[row][col][dep]   = totalVal + intZeroKrhoReal.tes[i][j][k] ;  break;
                        case InterMember::TMS: 
                            numIntSemi.tms[row][col][dep]   = totalVal + intZeroKrhoReal.tms[i][j][k] ;  break;
                        case InterMember::TER:  
                            numIntSemi.ter[row][col][dep]   = totalVal + intZeroKrhoReal.ter[i][j][k] ;  break;
                        case InterMember::TMR:  
                            numIntSemi.tmr[row][col][dep]   = totalVal + intZeroKrhoReal.tmr[i][j][k] ;  break;
                        case InterMember::TERZ1: 
                            numIntSemi.terz1[row][col][dep] = totalVal + intZeroKrhoReal.terz1[i][j][k]; break;
                        case InterMember::TMRZ1: 
                            numIntSemi.tmrz1[row][col][dep] = totalVal + intZeroKrhoReal.tmrz1[i][j][k]; break;
                        case InterMember::TERZ2: 
                            numIntSemi.terz2[row][col][dep] = totalVal + intZeroKrhoReal.terz2[i][j][k]; break;
                        case InterMember::TMRZ2: 
                            numIntSemi.tmrz2[row][col][dep] = totalVal + intZeroKrhoReal.tmrz2[i][j][k]; break;
                    }
                }
            }
        }
        
        // Imaginary axis integration path
        g_memoInter.clear();
        g_Data3D = dataImag3D;
        g_path = "imag";
        if (g_memoInter.empty()) g_memoInter.reserve(g_Data3D.indVecTriplet.size()); // Before integration loop (expected size)
        if (!g_Data3D.indVecTriplet.empty()) {
            IntegrandInter intZeroKrhoImag = numIntSemi.flatten(g_Data3D.indVecTriplet);
            // Loop over the members of the IntegrandInter structure
            for (InterMember mem : { InterMember::TE,    InterMember::TM,    InterMember::TEZ1,  InterMember::TMZ1,
                                     InterMember::TEZ2,  InterMember::TMZ2,  InterMember::TEZZ,  InterMember::TMZZ, 
                                     InterMember::TES,   InterMember::TMS,   InterMember::TER,   InterMember::TMR, 
                                     InterMember::TERZ1, InterMember::TMRZ1, InterMember::TERZ2, InterMember::TMRZ2 }) 
            {
                g_currentMemberInter = mem;
                for (size_t i = 0; i < g_Data3D.indVecTriplet.size(); i++) {
                    g_currentI = i; size_t j = 0; g_currentJ = j; size_t k = 0; g_currentK = k;
                    double realPart = std::real(gaussKronrod(pickRealMember, 1.0, 0.0, 5, 100, 1.0e-8, 1.0e-8));
                    double imagPart = std::real(gaussKronrod(pickImagMember, 1.0, 0.0, 5, 100, 1.0e-8, 1.0e-8));
                    dcmplx totalVal(realPart, imagPart);
                    int row = std::get<0>(g_Data3D.indVecTriplet[i]);
                    int col = std::get<1>(g_Data3D.indVecTriplet[i]);
                    int dep = std::get<2>(g_Data3D.indVecTriplet[i]);
                    switch (g_currentMemberInter)
                    {
                        case InterMember::TE:
                            numIntSemi.te[row][col][dep]    = totalVal + intZeroKrhoImag.te[i][j][k]  ;  break;
                        case InterMember::TM:
                            numIntSemi.tm[row][col][dep]    = totalVal + intZeroKrhoImag.tm[i][j][k]  ;  break;
                        case InterMember::TEZ1:
                            numIntSemi.tez1[row][col][dep]  = totalVal + intZeroKrhoImag.tez1[i][j][k] ; break;
                        case InterMember::TMZ1:
                            numIntSemi.tmz1[row][col][dep]  = totalVal + intZeroKrhoImag.tmz1[i][j][k] ; break;
                        case InterMember::TEZ2:
                            numIntSemi.tez2[row][col][dep]  = totalVal + intZeroKrhoImag.tez2[i][j][k] ; break;
                        case InterMember::TMZ2:
                            numIntSemi.tmz2[row][col][dep]  = totalVal + intZeroKrhoImag.tmz2[i][j][k] ; break;
                        case InterMember::TEZZ: 
                            numIntSemi.tezz[row][col][dep]  = totalVal + intZeroKrhoImag.tezz[i][j][k];  break;
                        case InterMember::TMZZ:
                            numIntSemi.tmzz[row][col][dep]  = totalVal + intZeroKrhoImag.tmzz[i][j][k];  break;
                        case InterMember::TES: 
                            numIntSemi.tes[row][col][dep]   = totalVal + intZeroKrhoImag.tes[i][j][k] ;  break;
                        case InterMember::TMS: 
                            numIntSemi.tms[row][col][dep]   = totalVal + intZeroKrhoImag.tms[i][j][k] ;  break;
                        case InterMember::TER:  
                            numIntSemi.ter[row][col][dep]   = totalVal + intZeroKrhoImag.ter[i][j][k] ;  break;
                        case InterMember::TMR:  
                            numIntSemi.tmr[row][col][dep]   = totalVal + intZeroKrhoImag.tmr[i][j][k] ;  break;
                        case InterMember::TERZ1: 
                            numIntSemi.terz1[row][col][dep] = totalVal + intZeroKrhoImag.terz1[i][j][k]; break;
                        case InterMember::TMRZ1: 
                            numIntSemi.tmrz1[row][col][dep] = totalVal + intZeroKrhoImag.tmrz1[i][j][k]; break;
                        case InterMember::TERZ2: 
                            numIntSemi.terz2[row][col][dep] = totalVal + intZeroKrhoImag.terz2[i][j][k]; break;
                        case InterMember::TMRZ2: 
                            numIntSemi.tmrz2[row][col][dep] = totalVal + intZeroKrhoImag.tmrz2[i][j][k]; break;
                    }
                }
            }
        }
        
#ifdef DEBUG
        if (materialIndex == 2 && sourceIndex == 1) {
            write3DMatrixToFile(numIntSemi.te,    "te.txt");    write3DMatrixToFile(numIntSemi.tm,    "tm.txt");
            write3DMatrixToFile(numIntSemi.tez1,  "tez1.txt");  write3DMatrixToFile(numIntSemi.tmz1,  "tmz1.txt");
            write3DMatrixToFile(numIntSemi.tez2,  "tez2.txt");  write3DMatrixToFile(numIntSemi.tmz2,  "tmz2.txt");
            write3DMatrixToFile(numIntSemi.tezz,  "tezz.txt");  write3DMatrixToFile(numIntSemi.tmzz,  "tmzz.txt");
            write3DMatrixToFile(numIntSemi.tes,   "tes.txt");   write3DMatrixToFile(numIntSemi.tms,   "tms.txt");
            write3DMatrixToFile(numIntSemi.ter,   "ter.txt");   write3DMatrixToFile(numIntSemi.tmr,   "tmr.txt");
            write3DMatrixToFile(numIntSemi.terz1, "terz1.txt"); write3DMatrixToFile(numIntSemi.tmrz1, "tmrz1.txt");
            write3DMatrixToFile(numIntSemi.terz2, "terz2.txt"); write3DMatrixToFile(numIntSemi.tmrz2, "tmrz2.txt");
        }
#endif

        CoeffsInter coeffs;
        std::vector<std::vector<std::vector<double>>> x = r;
        std::vector<std::vector<std::vector<double>>> y = z1;
        std::vector<std::vector<std::vector<double>>> z = z2;

        std::vector<double> xLims = getRowAtColumnDepth(x, 0, 0);
        std::vector<double> yLims = getColumnAtRowDepth(y, 0, 0);
        std::vector<double> zLims = getDepthAtRowColumn(z, 0, 0);

        double xMin = *std::min_element(xLims.begin(), xLims.end());
        double xMax = *std::max_element(xLims.begin(), xLims.end());
        double yMin = *std::min_element(yLims.begin(), yLims.end());
        double yMax = *std::max_element(yLims.begin(), yLims.end());
        double zMin = *std::min_element(zLims.begin(), zLims.end());
        double zMax = *std::max_element(zLims.begin(), zLims.end());

        std::vector<double> margins = {std::abs(xMax - xMin), std::abs(yMax - yMin), std::abs(zMax - zMin)};
        double maxScale = *std::max_element(margins.begin(), margins.end());
        std::vector<double> scales = {maxScale / margins[0], maxScale / margins[1], maxScale / margins[2]};
        coeffs.scales = scales;
        
        coeffs.te    = interpolateInter(x, y, z, scales, numIntSemi.te);    coeffs.tm    = interpolateInter(x, y, z, scales, numIntSemi.tm);
        coeffs.tez1  = interpolateInter(x, y, z, scales, numIntSemi.tez1);  coeffs.tmz1  = interpolateInter(x, y, z, scales, numIntSemi.tmz1);
        coeffs.tez2  = interpolateInter(x, y, z, scales, numIntSemi.tez2);  coeffs.tmz2  = interpolateInter(x, y, z, scales, numIntSemi.tmz2);
        coeffs.tezz  = interpolateInter(x, y, z, scales, numIntSemi.tezz);  coeffs.tmzz  = interpolateInter(x, y, z, scales, numIntSemi.tmzz);
        coeffs.tes   = interpolateInter(x, y, z, scales, numIntSemi.tes);   coeffs.tms   = interpolateInter(x, y, z, scales, numIntSemi.tms);
        coeffs.ter   = interpolateInter(x, y, z, scales, numIntSemi.ter);   coeffs.tmr   = interpolateInter(x, y, z, scales, numIntSemi.tmr);
        coeffs.terz1 = interpolateInter(x, y, z, scales, numIntSemi.terz1); coeffs.tmrz1 = interpolateInter(x, y, z, scales, numIntSemi.tmrz1);
        coeffs.terz2 = interpolateInter(x, y, z, scales, numIntSemi.terz2); coeffs.tmrz2 = interpolateInter(x, y, z, scales, numIntSemi.tmrz2);
      
        return coeffs;
    }
}

// Interface Fresnel reflection and transmission coefficients (TE or TM).
std::vector<dcmplx> SommerfeldIntegrator::FresnelCoeff(int layerIndexI, int layerIndexJ, 
                                                       dcmplx krho, const std::string& waves) const {
  
  // Validate layer indices
  if (layerIndexI < 0 || layerIndexI >= epsilon.size() || 
      layerIndexJ < 0 || layerIndexJ >= epsilon.size() || 
      std::abs(layerIndexI - layerIndexJ) > 1){
    throw std::out_of_range("Invalid layer indeces for Fresnel reflection coefficient calculation.");
  }

  // Check if the wave type is valid (TE or TM)
  if (waves != "TE" && waves != "TM"){
    throw std::invalid_argument("Invalid wave type. Use 'TE' or 'TM'.");
  }

  // Get the permittivity (TM) and permeability (TE) of the two layers
  dcmplx pI(0,0);
  dcmplx pJ(0,0);
  if (waves == "TM"){
    pI = epsilon[layerIndexI];
    pJ = epsilon[layerIndexJ];
  }
  else if (waves == "TE"){
    pI = mu[layerIndexI];
    pJ = mu[layerIndexJ];
  }
  
  // z component of the wavevector
  dcmplx kIz = k_L[layerIndexI];
  dcmplx kJz = k_L[layerIndexJ];
  kIz = kIz * kIz - krho * krho; kIz = csqrt(kIz);
  kJz = kJz * kJz - krho * krho; kJz = csqrt(kJz);

  // Fresnel reflection coefficient
  dcmplx fresnelReflCoeff;
  fresnelReflCoeff = (pJ * kIz - pI * kJz) /
                     (pJ * kIz + pI * kJz);

  // Fresnel transmission coefficient
  dcmplx fresnelTransCoeff;
  fresnelTransCoeff = (2.0 * pJ * kIz) / ( pJ * kIz + pI * kJz );

  std::vector<dcmplx> fresnelCoeffs = {fresnelReflCoeff, fresnelTransCoeff};
  
  return fresnelCoeffs;
}

// Write 2D matrix to .txt file for debugging purposes.
void SommerfeldIntegrator::writeMatrixToFile(const std::vector<std::vector<dcmplx>>& matrix, const std::string& filename) {
    std::ofstream file(filename);
    if (!file) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    for (size_t i = 0; i < matrix.size(); i++) {
        for (size_t j = 0; j < matrix[i].size(); j++) {
            file << i << " " << j << " " 
                 << real(matrix[i][j]) << " " << imag(matrix[i][j]) << "\n";
        }
    }
    file.close();
}

// Write 3D matrix to .txt file for debugging purposes.
void SommerfeldIntegrator::write3DMatrixToFile(const std::vector<std::vector<std::vector<dcmplx>>>& matrix, const std::string& filename) {
    std::ofstream file(filename);
    if (!file) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    for (size_t i = 0; i < matrix.size(); i++) {
        for (size_t j = 0; j < matrix[i].size(); j++) {
            for (size_t k = 0; k < matrix[i][j].size(); k++) {
                file << i << " " << j << " " << k << " "
                     << real(matrix[i][j][k]) << " " << imag(matrix[i][j][k]) << "\n";
            }
        }
    }
    file.close();
}
