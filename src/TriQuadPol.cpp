/**
 * @file TriQuadPol.cpp
 * @brief Polar-quadrature assembly about an observation point on a triangle.
 */

#include "TriQuadPol.h"

/**
 * @details Stores the inputs and calls adaptPolar() to generate nodes/weights.
 */
TriQuadPol::TriQuadPol(Triangle* inTp, rvec inPos, int nQuadrature,
                       double(*xQuadrature)[2], double* wQuadrature)
  : Tp(inTp), pos(inPos), nQuadrature(nQuadrature), 
    xQuadrature(xQuadrature), wQuadrature(wQuadrature) { adaptPolar(); }

/**
 * @details Main driver: project to plane, split angular sectors, and assemble nodes/weights.
 *
 * Clears previous results, then:
 *  1) transformInPlane()
 *  2) Decide inside/outside and compute sector data via adaptIn()/adaptOut()
 *  3) assembleFromPolar()
 */
void TriQuadPol::adaptPolar()
{    
    quadPoints.clear();     // clear previous results
    quadWeights.clear();    // clear previous results
    
    transformInPlane();     // project from 3D to 2D triangle plane

    // Center of the triangle in plane
    rvec TpCenter(0.0, 0.0, 0.0);
    for (int i = 0; i < 3; ++i) {
        TpCenter += vertsInPlane[i] / 3.0;
    }

    // Angle between impact positions and boundary centroids
    rvec distVec = TpCenter - posInPlane;
    double vertsTheta0 = std::atan2(distVec[1], distVec[0]);

    // Expand vertices in rotated system
    std::vector<double> vertX(3), vertY(3);
    for (int i = 0; i < 3; ++i) {
        vertX[i] = (vertsInPlane[i][0] - posInPlane[0]) * cos(vertsTheta0) +
                   (vertsInPlane[i][1] - posInPlane[1]) * sin(vertsTheta0);
        vertY[i] = (vertsInPlane[i][1] - posInPlane[1]) * cos(vertsTheta0) -
                   (vertsInPlane[i][0] - posInPlane[0]) * sin(vertsTheta0);
    }

    // Convert vertices to polar coordinates
    std::vector<double> vertsRho(3), vertsTheta(3);
    for (int i = 0; i < 3; ++i) {
        vertsRho[i] = std::hypot(vertX[i], vertY[i]);
        vertsTheta[i] = std::atan2(vertY[i], vertX[i]);
    }

    // Sort vertices by angle
    std::vector<int> idx(3);
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](int a, int b) { return vertsTheta[a] < vertsTheta[b]; });
    std::vector<double> vertsRhoSorted(3), vertsThetaSorted(3);
    for (int i = 0; i < 3; ++i) {
        vertsRhoSorted[i]   = vertsRho[idx[i]];
        vertsThetaSorted[i] = vertsTheta[idx[i]];
    }
    
    // Adapt integration points to triangle
    std::vector<double> intThetas, intWeights, minRadius, maxRadius;
    bool isInTriangle = pointInTriangle(posInPlane, vertsInPlane[0], vertsInPlane[1], vertsInPlane[2]);
    if (isInTriangle)    // if the observation point is within the triangle
    {
        // Line in polar coordinates: r(t) = d0 / cos(t - t0) -> Output: d0, t0
        std::pair<double, double> line01 = lineInPolar(vertsThetaSorted[0], vertsRhoSorted[0],
                                                       vertsThetaSorted[1], vertsRhoSorted[1]);
        std::pair<double, double> line12 = lineInPolar(vertsThetaSorted[1], vertsRhoSorted[1],
                                                       vertsThetaSorted[2], vertsRhoSorted[2]);
        std::pair<double, double> line20 = lineInPolar(vertsThetaSorted[2], vertsRhoSorted[2],
                                                       vertsThetaSorted[0], vertsRhoSorted[0]);

        // Adapt integration points to triangle
        std::vector<double> thetaTemp, weightTemp, minRadiusTemp, maxRadiusTemp;
        adaptIn(vertsTheta0, vertsThetaSorted[0], vertsThetaSorted[1], line01.first, line01.second,
                thetaTemp, weightTemp, minRadiusTemp, maxRadiusTemp);
        intThetas.insert(intThetas.end(), thetaTemp.begin(), thetaTemp.end());
        intWeights.insert(intWeights.end(), weightTemp.begin(), weightTemp.end());
        minRadius.insert(minRadius.end(), minRadiusTemp.begin(), minRadiusTemp.end());
        maxRadius.insert(maxRadius.end(), maxRadiusTemp.begin(), maxRadiusTemp.end());

        adaptIn(vertsTheta0, vertsThetaSorted[1], vertsThetaSorted[2], line12.first, line12.second,
                thetaTemp, weightTemp, minRadiusTemp, maxRadiusTemp);
        intThetas.insert(intThetas.end(), thetaTemp.begin(), thetaTemp.end());
        intWeights.insert(intWeights.end(), weightTemp.begin(), weightTemp.end());
        minRadius.insert(minRadius.end(), minRadiusTemp.begin(), minRadiusTemp.end());
        maxRadius.insert(maxRadius.end(), maxRadiusTemp.begin(), maxRadiusTemp.end());

        adaptIn(vertsTheta0, vertsThetaSorted[2], vertsThetaSorted[0], line20.first, line20.second,
                thetaTemp, weightTemp, minRadiusTemp, maxRadiusTemp);
        intThetas.insert(intThetas.end(), thetaTemp.begin(), thetaTemp.end());
        intWeights.insert(intWeights.end(), weightTemp.begin(), weightTemp.end());
        minRadius.insert(minRadius.end(), minRadiusTemp.begin(), minRadiusTemp.end());
        maxRadius.insert(maxRadius.end(), maxRadiusTemp.begin(), maxRadiusTemp.end());
    }
    else    // if the observation point is not within the circumcircle
    {
        // Line in polar coordinates: r(t) = d0 / cos(t - t0) -> Output: d0, t0
        std::pair<double, double> line01 = lineInPolar(vertsThetaSorted[0], vertsRhoSorted[0],
                                                       vertsThetaSorted[1], vertsRhoSorted[1]);
        std::pair<double, double> line02 = lineInPolar(vertsThetaSorted[0], vertsRhoSorted[0],
                                                       vertsThetaSorted[2], vertsRhoSorted[2]);
        std::pair<double, double> line12 = lineInPolar(vertsThetaSorted[1], vertsRhoSorted[1],
                                                       vertsThetaSorted[2], vertsRhoSorted[2]);

        // Adapt integration points to triangle
        std::vector<double> thetaTemp, weightTemp, minRadiusTemp, maxRadiusTemp;
        adaptOut(vertsTheta0, vertsThetaSorted[0], vertsThetaSorted[1],
                 line01.first, line01.second, line02.first, line02.second,
                 thetaTemp, weightTemp, minRadiusTemp, maxRadiusTemp);
        intThetas.insert(intThetas.end(), thetaTemp.begin(), thetaTemp.end());
        intWeights.insert(intWeights.end(), weightTemp.begin(), weightTemp.end());
        minRadius.insert(minRadius.end(), minRadiusTemp.begin(), minRadiusTemp.end());
        maxRadius.insert(maxRadius.end(), maxRadiusTemp.begin(), maxRadiusTemp.end());

        adaptOut(vertsTheta0, vertsThetaSorted[1], vertsThetaSorted[2],
                 line12.first, line12.second, line02.first, line02.second,
                 thetaTemp, weightTemp, minRadiusTemp, maxRadiusTemp);
        intThetas.insert(intThetas.end(), thetaTemp.begin(), thetaTemp.end());
        intWeights.insert(intWeights.end(), weightTemp.begin(), weightTemp.end());
        minRadius.insert(minRadius.end(), minRadiusTemp.begin(), minRadiusTemp.end());
        maxRadius.insert(maxRadius.end(), maxRadiusTemp.begin(), maxRadiusTemp.end());
    }
    // Assemble quadrature points
    assembleFromPolar(intThetas, intWeights, minRadius, maxRadius);
}

/**
 * @details Builds an orthonormal basis (e1 along edge 0->1, e2 = e3 x e1, e3 normal),
 * then expresses vertices and the observation point in that basis.
 */
void TriQuadPol::transformInPlane() {
    rvec rp[3];
    for (int i = 0; i < 3; ++i) rp[i] = Tp->Node(i);

    // unit vectors
    rvec e1 = rp[1] - rp[0]; e1 /= sqrt( dot(e1, e1) );
    rvec temp = rp[2] - rp[0];
    rvec e3 = cross(e1, temp); e3 /= sqrt( dot(e3, e3) );
    rvec e2 = cross(e3, e1);

    // vertices in triangle plane
    vertsInPlane.push_back( rvec( 0.0, 0.0, 0.0 ) ); // rp[0] in plane
    temp = rp[1] - rp[0];
    vertsInPlane.push_back( rvec( clean( dot(temp, e1) ), 
                                  clean( dot(temp, e2) ), 0.0 ) ); // rp[1] in plane
    temp = rp[2] - rp[0];
    vertsInPlane.push_back( rvec( clean( dot(temp, e1) ), 
                                  clean( dot(temp, e2) ), 0.0 ) ); // rp[2] in plane

    // Origin in triangle plane
    posInPlane = rvec( clean( dot(e1, pos - rp[0]) ), clean( dot(e2, pos - rp[0]) ), 0.0 );
}

/**
 * @details Return line parameters (θ0,d0) for the edge between two vertices in polar coordinates.
 *
 * To avoid degeneracy when θ1 == θ2, a tiny perturbation is applied.
 */
std::pair<double, double> TriQuadPol::lineInPolar(double theta1, double rho1, double theta2, double rho2) {
    if (theta1 == theta2) theta2 += 1e-10; // avoid division by zero
    double sinSign = std::sin(theta2 - theta1) > 0 ? 1.0 : -1.0;
    double R = std::sqrt(rho1 * rho1 + rho2 * rho2 - 2.0 * rho1 * rho2 * std::cos(theta2 - theta1));

    // line parameters
    double coeff0 = sinSign / R * rho1 * rho2 * std::sin(theta2 - theta1);
    double theta0 = std::atan2(-sinSign * (rho2 * std::cos(theta2) - rho1 * std::cos(theta1)),
                                sinSign * (rho2 * std::sin(theta2) - rho1 * std::sin(theta1)) );
    return std::make_pair(theta0, coeff0);
}

/**
 * @details Build a θ-bin and radial bound for the inside case (single edge).
 *
 * Fills @p intThetas/intWeights with Gauss nodes/weights on [θ_start, θ_end],
 * and sets r_max(θ) = d0 / cos(θ − θ0), r_min = 0.
 */
void TriQuadPol::adaptIn(double theta0, double thetaStart, double thetaEnd, double tLine, double dLine,
                         std::vector<double>& intThetas, std::vector<double>& intWeights,
                         std::vector<double>& minRadius, std::vector<double>& maxRadius)
{
    // clear previous results
    intThetas.clear(); intWeights.clear(); minRadius.clear(); maxRadius.clear();

    // Wrap-around for negative angle spans
    if (thetaEnd < thetaStart) thetaEnd += 2.0 * M_PI;

    // Allocate and compute
    intThetas.resize(nQuadrature); intWeights.resize(nQuadrature);
    minRadius.resize(nQuadrature, 0.0); maxRadius.resize(nQuadrature);

    for (int i = 0; i < nQuadrature; ++i) {
        intThetas[i] = thetaStart + (thetaEnd - thetaStart) * xQuadrature[i][0];
        maxRadius[i] = dLine / std::cos(intThetas[i] - tLine);
        intWeights[i] = (thetaEnd - thetaStart) * wQuadrature[i];
        intThetas[i] += theta0;
    }
}

/**
 * @details Build a θ-bin and two radial bounds for the outside case (two edges).
 *
 * Sets r_min/max(θ) from two edge lines; swaps if min > max numerically.
 */
void TriQuadPol::adaptOut(double theta0, double thetaStart, double thetaEnd,
                          double tLineFirst, double dLineFirst, double tLineSecond, double dLineSecond,
                          std::vector<double>& intThetas, std::vector<double>& intWeights,
                          std::vector<double>& minRadius, std::vector<double>& maxRadius)
{
    // clear previous results
    intThetas.clear(); intWeights.clear(); minRadius.clear(); maxRadius.clear();

    // Allocate and compute
    intThetas.resize(nQuadrature); intWeights.resize(nQuadrature);
    minRadius.resize(nQuadrature, 0.0); maxRadius.resize(nQuadrature);

    for (int i = 0; i < nQuadrature; ++i) {
        intThetas[i] = thetaStart + (thetaEnd - thetaStart) * xQuadrature[i][0];
        minRadius[i] = dLineFirst  / std::cos(intThetas[i] - tLineFirst);
        maxRadius[i] = dLineSecond / std::cos(intThetas[i] - tLineSecond);
        intWeights[i] = (thetaEnd - thetaStart) * wQuadrature[i];
        intThetas[i] += theta0;
        if (minRadius[i] > maxRadius[i]) std::swap(minRadius[i], maxRadius[i]);
    }
}

/**
 * @details Expand (θ,r) Gauss products and convert to barycentric coordinates.
 *
 * Each θ-bin is expanded by nQuadrature radial nodes. The Jacobian factor r
 * and the radial span (r_max - r_min) are included in the weights.
 */
void TriQuadPol::assembleFromPolar(const std::vector<double>& inThetas, 
                                   const std::vector<double>& inWeights,
                                   const std::vector<double>& inMinRadii,
                                   const std::vector<double>& inMaxRadii)
{
    // triangle vertices
    const rvec& vertex1 = vertsInPlane[0];
    const rvec& vertex2 = vertsInPlane[1];
    const rvec& vertex3 = vertsInPlane[2];

    for (int j = 0; j < inThetas.size(); ++j) {
        for (int i = 0; i < nQuadrature; ++i) {
            // expand integration points and weights along radial direction
            double r = inMinRadii[j] + xQuadrature[i][0] * (inMaxRadii[j] - inMinRadii[j]);
            double weight = inWeights[j] * wQuadrature[i] * (inMaxRadii[j] - inMinRadii[j]) * r;
            
            // convert integration points to Cartesian coordinates
            double x = posInPlane[0] + r * std::cos(inThetas[j]);
            double y = posInPlane[1] + r * std::sin(inThetas[j]);

            // convert from Cartesian to triangle coordinates
            double a = (vertex2[1] - vertex3[1]) * (vertex1[0] - vertex3[0]) + 
                       (vertex3[0] - vertex2[0]) * (vertex1[1] - vertex3[1]);
            double xt = ( (vertex2[1] - vertex3[1]) * (x - vertex3[0]) + 
                          (vertex3[0] - vertex2[0]) * (y - vertex3[1]) ) / a;
            double yt = ( (vertex3[1] - vertex1[1]) * (x - vertex3[0]) + 
                          (vertex1[0] - vertex3[0]) * (y - vertex3[1])) / a;

            rvec pt(3); pt[0] = xt; pt[1] = yt; pt[2] = 1.0 - xt - yt;
            quadPoints.push_back(pt);
            quadWeights.push_back(weight);
        }
    }
}