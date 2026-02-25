/**
 * @file SingSub.cpp
 * @brief Implementation of exact edge/surface integrals and K1–K4 kernels.
 */

#include "SingSub.h"
#include <cmath>
#include "GreenF.h"

/**
 * @def SURFTHRESH
 * @brief Numerical threshold for detecting points lying near an edge or plane.
 *
 * Used in geometric tests to treat points within @c 1e-8 of a surface
 * as being on the surface.
 */
#define SURFTHRESH 1e-8

/** @brief Default constructor. */
SingSub::SingSub() {}

/**
 * @details Initialize with triangle from RWG and observation point r.
 *
 * Computes:
 *  - unit edge-normal directions (tangent x RefNormalVector) for each edge,
 *  - signed distances t_i from r to the three supporting lines,
 *  - signed plane distance h = n(r − r0), with n = RefNormalVector(),
 *  - signed solid angle omega subtended by the triangle at r.
 *
 * The solid angle sign follows the sign of h.
 */
int SingSub::Init(RWGFun* rwgIn, rvec rIn) {
  r = rIn;

  // Check if triangle has changed, only then recompute edge normals
  if (rwgIn != rwgFun) {
    rwgFun = rwgIn;
    edgeNormals.clear();

    // Edge normals: e_i = (v_{i+1} - v_i) x n_ref, normalized
    edgeNormals.push_back(cross(
        (rvec)(rwgFun->TrianglePtr()->Node(1) - rwgFun->TrianglePtr()->Node(0)),
        rwgFun->TrianglePtr()->RefNormalVector()));
    edgeNormals.push_back(cross(
        (rvec)(rwgFun->TrianglePtr()->Node(2) - rwgFun->TrianglePtr()->Node(1)),
        rwgFun->TrianglePtr()->RefNormalVector()));
    edgeNormals.push_back(cross(
        (rvec)(rwgFun->TrianglePtr()->Node(0) - rwgFun->TrianglePtr()->Node(2)),
        rwgFun->TrianglePtr()->RefNormalVector()));

    // Normalize
    for (int i = 0; i < 3; i++)
      edgeNormals[i] /= sqrt(dot(edgeNormals[i], edgeNormals[i]));
  }

  // Calculate "t" factors
  t.clear();
  for (int i = 0; i < 3; i++) {
    t.push_back(
        dot(edgeNormals[i], (rvec)(r - rwgFun->TrianglePtr()->Node(i))));
  }

  // Calculate h and omega
  h = dot(rwgFun->TrianglePtr()->RefNormalVector(),
          (rvec)(r - rwgFun->TrianglePtr()->Node(0)));

  // Auxiliary unit vectors from r to vertices
  rvec a1, a2, a3;
  double fNorm;

  a1 = rwgFun->TrianglePtr()->Node(0) - r;
  fNorm = sqrt(dot(a1, a1));
  if (fNorm > 0) a1 /= fNorm;

  a2 = rwgFun->TrianglePtr()->Node(1) - r;
  fNorm = sqrt(dot(a2, a2));
  if (fNorm > 0) a2 /= fNorm;

  a3 = rwgFun->TrianglePtr()->Node(2) - r;
  fNorm = sqrt(dot(a3, a3));
  if (fNorm > 0) a3 /= fNorm;

  double x = 1. + dot(a1, a2) + dot(a1, a3) + dot(a2, a3);
  double y = fabs(dot(a1, cross(a2, a3)));

  if (h > 0.)
    omega = 2. * fabs(atan2(y, x));
  else
    omega = -2. * fabs(atan2(y, x));
  return 0;
}

/**
 * @details Preferred initialization: precompute ILar/ISar for the triangle at r.
 *
 * Sets:
 *  - h and omega as in Init(),
 *  - eNar[i] (unit edge-normals) and tar[i] = −eNar[i]·(Node(i)−r),
 *  - ILar[i][k] for q = {-1,1,3} (k = 0,1,2),
 *  - ISar[k]    for q = {-3,-1,1,3} (k = 0..3).
 */
int SingSub::SetTriangleRWG(RWGFun* rwgIn, rvec rIn) {
  r = rIn;
  rwgFun = rwgIn;
  tPtr = rwgFun->TrianglePtr();

  // Signed distance to plane
  h = dot(tPtr->RefNormalVector(), (rvec)(r - tPtr->Node(0)));

  // Build unit edge tangents and corresponding edge normals
  double fNorm;
  rvec edge[3], a[3], anorm[3];
  for (int i = 0; i < 3; ++i) {
    int ii = (i + 1) % 3;
    edge[i] = tPtr->Node(ii) - tPtr->Node(i);
    edge[i] /= sqrt(dot(edge[i], edge[i]));
    eNar[i] = cross(edge[i], tPtr->RefNormalVector());
    eNar[i] /= sqrt(dot(eNar[i], eNar[i]));

    a[i] = tPtr->Node(i) - r;
    fNorm = sqrt(dot(a[i], a[i]));
    if (fNorm > 0)
      anorm[i] = a[i] / fNorm;
    else
      anorm[i] = a[i];

    tar[i] = -dot(eNar[i], a[i]);
  }

  // Signed solid angle using unit vectors anorm[i]
  double x = 1. + dot(anorm[0], anorm[1]) + dot(anorm[1], anorm[2]) +
             dot(anorm[2], anorm[0]);
  double y = fabs(dot(anorm[0], cross(anorm[1], anorm[2])));

  if (h > 0.) {
    omega = 2. * fabs(atan2(y, x));
  } else {
    omega = -2. * fabs(atan2(y, x));
  }

  // Precompute line integrals on each edge for q = {-1,1,3}
  for (int i = 0; i < 3; ++i) {
    int ii = (i + 1) % 3;
    double sp(dot(edge[i], a[ii]));
    double sm(dot(edge[i], a[i]));

    double Rp(sqrt(dot(a[ii], a[ii])));
    double Rm(sqrt(dot(a[i], a[i])));
    double R0(Rp * Rp - sp * sp);

    if (R0 < 0.) R0 = 0.;

    R0 = sqrt(R0);

    // Log term with endpoint handling (Rp==0 or Rm==0)
    double fLog;
    if (Rp == 0 || Rm == 0)
      fLog = 0.;
    else if (Rm + sm > Rp - sp)
      fLog = log((Rp + sp) / (Rm + sm));
    else
      fLog = log((Rm - sm) / (Rp - sp));

    if (fabs(R0) < SURFTHRESH) {
      ILar[i][0] = 0;
      ILar[i][1] = (sp * Rp - sm * Rm) / 2.;
      ILar[i][2] = (sp * Rp * Rp * Rp - sm * Rm * Rm * Rm) / 4.;
    } else {
      ILar[i][0] = fLog;
      ILar[i][1] = (R0 * R0 * fLog + sp * Rp - sm * Rm) / 2.;
      ILar[i][2] = R0 * R0 * 3. / 8. * (R0 * R0 * fLog + sp * Rp - sm * Rm) +
                   (sp * Rp * Rp * Rp - sm * Rm * Rm * Rm) / 4.;
    }
  }

  // Precompute surface integrals ISar[q]
  ISar[0] = omega / h;
  ISar[1] = omega * h;
  ISar[2] = omega * h * h * h;
  ISar[3] = omega * h * h * h * h * h;
  for (int i = 0; i < 3; ++i) {
    if (fabs(tar[i]) > SURFTHRESH) {
      ISar[1] += tar[i] * ILar[i][0];
      ISar[2] += tar[i] * (h * h * ILar[i][0] + ILar[i][1]);
      ISar[3] += tar[i] *
                 (h * h * h * h * ILar[i][0] + h * h * ILar[i][1] + ILar[i][2]);
    }
  }
  ISar[1] /= -1;
  ISar[2] /= -3;
  ISar[3] /= -5;

  return 0;
}

// Change only the RWG/triangle while keeping r fixed.
int SingSub::ChangeTriangleRWG(RWGFun* rwgIn) {
  rwgFun = rwgIn;
  return 0;
}

// Scalar K1 term
double SingSub::K1(int q) { return ISar[(q + 3) / 2] * rwgFun->RWGPreFactor(); }

// Vector K2 term
rvec SingSub::K2(int q) {
  rvec rho(r - h * tPtr->RefNormalVector());
  rvec result((rho - rwgFun->FreeVertex()) * ISar[(q + 3) / 2]);
  for (int i = 0; i < 3; i++) {
    result += eNar[i] * ILar[i][(q + 3) / 2] / (q + 2.);
  }
  result *= rwgFun->RWGPreFactor() / 2.;
  return result;
}

// Vector K3 term
rvec SingSub::K3(int q) {
  rvec result(0., 0., 0.);
  for (int i = 0; i < 3; i++) result += eNar[i] * ILar[i][(q + 1) / 2];
  if (q == -1)
    result += tPtr->RefNormalVector() * omega;
  else
    result -= h * q * tPtr->RefNormalVector() * ISar[(q + 1) / 2];
  return result;
}

// Normal-only part of K3(q)
rvec SingSub::K3n(int q) {
  if (q == -1)
    return tPtr->RefNormalVector() * omega;
  else
    return -h * q * tPtr->RefNormalVector() * ISar[(q + 1) / 2];
}

// Vector K4 term
rvec SingSub::K4(int q) {
  return cross((rvec)(r - rwgFun->FreeVertex()), K3(q)) *
         rwgFun->RWGPreFactor() / (-2.);
}

/**
 * @details Exact line integral of @f$ R^q @f$ along edge triEdgeIndex.
 *
 * Handles endpoints where \p Rp==0 or \p Rm==0 in a way suitable for singularity
 * subtraction: the logarithmic divergence of IL(-1) at the vertex is cancelled
 * by its prefactor in the combined expressions; for q>−1, the log term is
 * multiplied by R0^2 which vanishes at those endpoints.
 */
double SingSub::IL(int q, int triEdgeIndex) {
  rvec p1(rwgFun->TrianglePtr()->Node(triEdgeIndex % 3));
  rvec p2(rwgFun->TrianglePtr()->Node((triEdgeIndex + 1) % 3));
  rvec vTemp(p2 - p1);
  rvec s(vTemp / sqrt(dot(vTemp, vTemp)));
  double sp(dot(s, (rvec)(p2 - r)));
  double sm(dot(s, (rvec)(p1 - r)));

  vTemp = r - p2;
  double Rp(sqrt(dot(vTemp, vTemp)));
  vTemp = r - p1;
  double Rm(sqrt(dot(vTemp, vTemp)));
  double R0(Rp * Rp - sp * sp);

  if (R0 < 0.) R0 = 0.;

  R0 = sqrt(R0);

  double fLog;
  // Note: this is not strictly correct! For r on either p1 or p2,
  // IL(-1) will diverge logarithmically. However, in the singularity
  // subtraction technique only t_i * IL(-1) appears, the result will
  // be equal to zero and so this is acceptable. Also, for q > -1,
  // fLog is multiplied by R0^2, which is equal to zero for these
  // cases, so again, this is acceptable:
  if (Rp == 0 || Rm == 0)
    fLog = 0.;
  else if (Rm + sm > Rp - sp)
    fLog = log((Rp + sp) / (Rm + sm));
  else
    fLog = log((Rm - sm) / (Rp - sp));

  switch (q) {
    case -1:
      return fLog;
    case 1:
      if (fabs(R0) < SURFTHRESH)
        return (sp * Rp - sm * Rm) / 2.;
      else
        return (R0 * R0 * fLog + sp * Rp - sm * Rm) / 2.;
    case 3:
      if (fabs(R0) < SURFTHRESH)
        return (sp * Rp * Rp * Rp - sm * Rm * Rm * Rm) / 4.;
      else
        return R0 * R0 * 3. / 8. * (R0 * R0 * fLog + sp * Rp - sm * Rm) +
               (sp * Rp * Rp * Rp - sm * Rm * Rm * Rm) / 4.;
  }
  return 0.;
}

// Exact surface integral of @f$ R^q @f$ over the triangle.
double SingSub::IS(int q) {
  double result(0);
  switch (q) {
    case -3:
      result = omega / h;
      break;
    case -1:
      result = h * omega;
      for (int i = 0; i < 3; i++) {
        if (fabs(t[i]) > SURFTHRESH) result += t[i] * IL(-1, i);
      }
      result *= -1.;
      break;
    case 1:
      result = h * h * h * omega;
      for (int i = 0; i < 3; i++)
        if (fabs(t[i]) > SURFTHRESH)
          result += t[i] * (h * h * IL(-1, i) + IL(1, i));
      result /= -3.;
      break;
    case 3:
      result = h * h * h * h * h * omega;
      for (int i = 0; i < 3; i++)
        if (fabs(t[i]) > SURFTHRESH)
          result +=
              t[i] * (h * h * h * h * IL(-1, i) + h * h * IL(1, i) + IL(3, i));
      result /= -5.;
  }
  return result;
}

