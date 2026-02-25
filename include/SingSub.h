/**
 * @file SingSub.h
 * @brief Exact subtractions for singular/near-singular triangle integrals.
 *
 * The SingSub class provides *closed-form* line/surface integrals of powers of
 * the distance R between an observation point r and points on a triangle.
 * These are used to compute the analytic parts in singularity subtraction,
 * and to form the vector/tensor K1–K4 terms. This implementation follows
 * Hänninen, I., Taskinen, M., & Sarvas, J. (2006)
 */

// Only include this once during compiling
#ifndef SINGSUB_H
#define SINGSUB_H

#include <vector>
#include "RWGFun.h"
#include "Triangle.h"
#include "globals.h"

/// \ingroup quadrature
/**
 * @class SingSub
 * @brief Analytic helpers for singularity subtraction on a single triangle.
 *
 * Given a triangle (from a RWG function) and an observation point r, this
 * class computes line and surface integrals of @f$ R^q (q \in \{-3,-1,1,3\}) @f$, 
 * and the derived vector expressions \p K1-K4 used in PMCHWT kernels.
 *
 * Typical use:
 *   - Call SetTriangleRWG() once for a new triangle/point pair.
 *   - Use \p K1-K4 or \p IL/IS repeatedly without recomputing geometry.
 *   - If only the RWG/triangle changes but r stays the same, call
 *     ChangeTriangleRWG() and reuse the cached r-dependent terms.
 */
class SingSub {
 private:
  // --- Cached state (from RWG/triangle and observation point) ---
  RWGFun* rwgFun;                 ///< RWG function defining the triangle
  rvec r;                         ///< Observation point
  std::vector<rvec> edgeNormals;  ///< Unit normals of triangle's edges
  std::vector<double> t;          ///< Signed distances to the three supporting lines
  double h;                       ///< Signed distance from r to the triangle plane
  double omega;                   ///< Signed solid angle subtended by the triangle at r

  // --- Precomputed arrays for "better" initialization path ---
  Triangle* tPtr;     ///< Convenience pointer to the triangle
  rvec eNar[3];       ///< Edge-normal directions (unit)
  double tar[3];      ///< -dot(eNar[i], a[i]) with a[i] = Node(i) − r
  double ILar[3][3];  ///< Edge integrals (edge index 0...2; columns for q = -1, 1, 3 -> indices 0,1,2)
  double ISar[4];     ///< Surface integrals for q = -3, -1, 1, 3 -> indices 0..3

 public:
  /** @brief Default constructor. */
  SingSub();

  /**
   * @brief Initialize with triangle from RWG and observation point r.
   *
   * Computes edge normals, distances t_i, signed plane distance h, and
   * signed solid angle omega. Use this when you want a single call that
   * refreshes all geometric cache for both triangle and r.
   *
   * @param rwgIn RWG function whose triangle is integrated over.
   * @param rIn   Observation point.
   * @return 0 on success.
   */
  int Init(RWGFun* rwgIn, rvec rIn);

  /**
   * @brief Preferred initialization that fills all cached arrays (ILar/ISar).
   *
   * Precomputes the line integrals along edges and the surface integrals
   * for @f$q \in \{-3,-1,1,3\}@f$, saving them in \p ILar/ISar; subsequent calls to
   * \p K1-K4/IL/IS reuse these values.
   *
   * @param rwgIn RWG function defining the triangle.
   * @param rIn   Observation point.
   * @return 0 on success.
   */
  int SetTriangleRWG(RWGFun* rwgIn, rvec rIn);

  /**
   * @brief Update only the RWG/triangle; keep the same observation point r.
   * @param rwgIn New RWG function.
   * @return 0 on success.
   */
  int ChangeTriangleRWG(RWGFun* rwgIn);

  /**
   * @brief Exact line integral \f$\int_{edge} R^q\,ds\f$ along a triangle edge.
   *
   * @param q            Power of R (allowed: -1, 1, 3).
   * @param triEdgeIndex Edge index (0,1,2) in the triangle's local ordering.
   * @return Value of the integral.
   */
  double IL(int q, int triEdgeIndex);
  
  /**
   * @brief Exact surface integral \f$\iint_{triangle} R^q\,dS\f$.
   * @param q Power of R (allowed: -3, -1, 1, 3).
   * @return Value of the integral.
   */
  double IS(int q);

  /**
   * @brief K1(q) scalar term.
   * @param q One of {-3, -1, 1, 3}.
   * @return Scalar K1 value.
   */
  double K1(int q);

  /**
   * @brief K2(q) vector term.
   * @param q One of {-3, -1, 1, 3}.
   * @return Vector K2 value.
   */
  rvec K2(int q);
  
  /**
   * @brief K3(q) vector term.
   * @param q One of {-3, -1, 1, 3}.
   * @return Vector K3 value.
   */
  rvec K3(int q);

  /**
   * @brief K3n(q) = normal-only contribution of K3(q).
   * @param q One of {-3, -1, 1, 3}.
   * @return Vector aligned with RefNormalVector().
   */
  rvec K3n(int q);

  /**
   * @brief K4(q) vector term.
   * @param q One of {-3, -1, 1, 3}.
   * @return Vector K4 value.
   */
  rvec K4(int q);
};

#endif
