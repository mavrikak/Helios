/**
 * @file TriQuadPol.h
 * @brief Polar quadrature on a single triangle around an observation point.
 *
 * Builds angular–radial integration segments in the triangle's plane,
 * then expands them via 1D Gauss rules to produce barycentric nodes and
 * weights suitable for integrating singular/near-singular kernels.
 */

#ifndef TRIQUADPOL_H
#define TRIQUADPOL_H

#include "Triangle.h"
#include "RWGFun.h"
#include "quadrature.h"
#include <algorithm>
#include <numeric>

/// \ingroup quadrature
/**
 * @class TriQuadPol
 * @brief Polar @f$(\rho, \theta)@f$ quadrature generator on a triangular facet.
 *
 * Workflow:
 *  - Construct with a triangle, an observation point, and a 1D Gauss rule
 *     (nodes/weights) of size @p nQuadrature for both @f$\rho@f$ and @f$\theta@f$ 
 *     expansions.
 *  - The constructor calls adaptPolar() which:
 *     + Projects geometry to the triangle plane (transformInPlane()).
 *     + Determines angular sectors and their radial bounds
 *       (adaptIn() if the point lies inside the facet's projection,
 *        else adaptOut()).
 *     + Assembles barycentric quadrature points/weights (assembleFromPolar()).
 *  - Retrieve integration nodes/weights via getQuadPoints()/#getQuadWeights().
 */
class TriQuadPol {
private:
    Triangle* Tp;                       ///< Triangle being integrated over
    rvec pos;                           ///< 3D observation point
    std::vector<rvec> vertsInPlane;     ///< Triangle vertices in the local 2D plane
    rvec posInPlane;                    ///< Observation point in the local 2D plane
    int nQuadrature;                    ///< Size of the 1D Gauss rule
    double(*xQuadrature)[2];            ///< 1D Gauss nodes (use [i][0]); [i][1] kept for compatibility
    double* wQuadrature;                ///< 1D Gauss weights
    std::vector<rvec> quadPoints;       ///< Output barycentric nodes @f$ (\lambda_1, \lambda_2, \lambda_3) @f$
    std::vector<double> quadWeights;    ///< Output weights (include Jacobians)

public:
    /**
     * @brief Construct and immediately adapt the polar quadrature.
     *
     * @param inTp         Triangle pointer.
     * @param inPos        Observation point in 3D.
     * @param nQuadrature  Number of nodes for the 1D Gauss rule (@f$\rho@f$ and @f$\theta@f$).
     * @param xQuadrature  Array of Gauss nodes in [0,1] (size nQuadrature x 2; uses [i][0]).
     * @param wQuadrature  Array of corresponding Gauss weights (size nQuadrature).
     */
    TriQuadPol(Triangle* inTp, rvec inPos, int nQuadrature,
               double(*xQuadrature)[2], double* wQuadrature);

    /**
     * @brief Build angular sectors and radial bounds, then assemble nodes/weights.
     *
     * Calls, in order: transformInPlane() -> adaptIn()/#adaptOut() -> assembleFromPolar().
     * Safe to call multiple times (clears previous results).
     */
    void adaptPolar();

    /**
     * @brief Build an orthonormal basis on the triangle plane and project (vertices, point).
     *
     * Sets @p vertsInPlane and @p posInPlane in a local (e1,e2) frame with vertex 0 at (0,0).
     */
    void transformInPlane();

    /**
     * @brief Return true if a, b, c all have the same sign (non-strict).
     * @param a,b,c Values to compare.
     * @return True if all have the same sign.
     */
    bool sameSign(double a, double b, double c) {
        return (a >= 0.0 && b >= 0.0 && c >= 0.0) ||
               (a <= 0.0 && b <= 0.0 && c <= 0.0);
    }

    /**
     * @brief 2D triangle membership test (local plane coordinates).
     * @param P Point to test.
     * @param A,B,C Triangle vertices.
     * @return true if P lies in triangle ABC (edge-inclusive).
     */
    bool pointInTriangle(const rvec& P, const rvec& A, const rvec& B, const rvec& C) {
        auto outerProd = [](const rvec& p1, const rvec& p2, const rvec& p3) {
            return (p1[0] - p3[0]) * (p2[1] - p3[1]) - (p2[0] - p3[0]) * (p1[1] - p3[1]);
        };

        double d1 = outerProd(P, A, B);
        double d2 = outerProd(P, B, C);
        double d3 = outerProd(P, C, A);

        return sameSign(d1, d2, d3);
    }

    /**
     * @brief Zero out tiny values for numeric stability.
     * @param x Value to clean.
     * @param eps Threshold.
     * @return Zero if |x| < eps, else x.
     */
    double clean(double x, double eps = 1e-10) {
        return std::abs(x) < eps ? 0.0 : x;
    }

    /**
     * @brief Straight-line segment between two vertices written in polar form about posInPlane.
     *
     * For a line through @f$(\theta_1,\rho_1)@f$ and @f$(\theta_2,\rho_2)@f$ in polar coords 
     * around posInPlane, returns @f$(\theta_0,\rho_0)@f$ such that 
     * @f$ r(\theta) = d_0 / cos(\theta − \theta_0)@f$ on that segment.
     *
     * @param theta1, rho1 Endpoint 1 in polar coords.
     * @param theta2, rho2 Endpoint 2 in polar coords.
     * @return Pair @f$\{\theta_0, d_0\}@f$.
     */
    std::pair<double, double> lineInPolar(double theta1, double rho1, double theta2, double rho2);

    /**
     * @brief Build angular bins and radial bounds when the point is inside the triangle.
     *
     * Uses a single boundary line per angular sector.
     *
     * @param theta0 Global rotation added to all output angles.
     * @param thetaStart, thetaEnd Sector angular limits.
     * @param tLine,dLine Line parameters from lineInPolar() @f$\{\theta_0, d_0\}@f$.
     * @param intThetas,intWeights Output Gauss nodes/weights in @f$\theta@f$.
     * @param minRadius,maxRadius  Output inner/outer radii for each @f$\theta@f$ bin.
     */
    void adaptIn(double theta0, double thetaStart, double thetaEnd, double tLine, double dLine,
                 std::vector<double>& intThetas, std::vector<double>& intWeights,
                 std::vector<double>& minRadius, std::vector<double>& maxRadius);

    /**
     * @brief Build angular bins/radial bounds when the point is outside the triangle.
     *
     * Uses two boundary lines per angular sector (enter/exit).
     *
     * @param theta0 Global rotation added to all output angles.
     * @param thetaStart, thetaEnd Sector angular limits.
     * @param tLineFirst,dLineFirst  First line parameters @f$\{\theta_0, d_0\}@f$.
     * @param tLineSecond,dLineSecond Second line parameters @f$\{\theta_0, d_0\}@f$.
     * @param intThetas,intWeights Output Gauss nodes/weights in @f$\theta@f$.
     * @param minRadius,maxRadius  Output inner/outer radii for each @f$\theta@f$ bin.
     */
    void adaptOut(double theta0, double thetaStart, double thetaEnd,
                  double tLineFirst, double dLineFirst, double tLineSecond, double dLineSecond,
                  std::vector<double>& intThetas, std::vector<double>& intWeights,
                  std::vector<double>& minRadius, std::vector<double>& maxRadius);

    /**
     * @brief Expand @f$\theta@f$-bins into (@f$\theta, r@f$) Gauss products and convert to barycentric nodes.
     *
     * @param inThetas   Vector of @f$\theta@f$ nodes (bin centers).
     * @param inWeights  Vector of @f$\theta@f$ weights.
     * @param inMinRadii Vector of @f$r_{min}@f$ per @f$\theta@f$ node.
     * @param inMaxRadii Vector of @f$r_{max}@f$ per @f$\theta@f$ node.
     */
    void assembleFromPolar(const std::vector<double>& inThetas, const std::vector<double>& inWeights,
                           const std::vector<double>& inMinRadii, const std::vector<double>& inMaxRadii);

    /** 
     * @brief Quadrature nodes (barycentric triples).
     * @return Vector of barycentric triples.
     */
    std::vector<rvec> getQuadPoints() { return quadPoints; }

    /**
     * @brief Quadrature weights (include all Jacobian factors and r).
     * @return Vector of weights.
     */
    std::vector<double> getQuadWeights() { return quadWeights; }
};

#endif
