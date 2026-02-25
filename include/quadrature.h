/**
 * @file quadrature.h
 * @brief Declarations of Dunavant (triangle) and Gauss–Legendre (1D) quadrature rules.
 *
 * Provides compile-time tables of nodes (barycentric / 1D) and weights for:
 *   - Dunavant Orders 1 (1 pt), 5 (7 pts), and 17 (61 pts) on triangles.
 *     + Nodes xN are given as barycentric triples @f$ (\lambda_1, \lambda_2, \lambda_3) @f$, row-major.
 *     + Weights wN integrate over the reference triangle with area 1.
 *   - Gauss–Legendre Orders 1, 3, 4, 8, and 17 on [0,1] (mapped from [-1,1]).
 *     + Nodes xN are given as @f$ (\xi, 1-\xi) @f$ pairs for convenience in tensor products.
 *     + Weights wN integrate over the unit interval [0,1].
 */

#ifndef QUADRATURE_H
#define QUADRATURE_H

/// \ingroup quadrature
namespace dunavant {
    /// Order-1 (1-point) Dunavant rule on a triangle: centroid only.
    extern int N1;
    extern double x1[1][3];   ///< Barycentric nodes, size N1x3
    extern double w1[1];      ///< Weights, size N1

    /// Order-5 (7-point) Dunavant rule on a triangle.
    extern int N5;
    extern double x5[7][3];   ///< Barycentric nodes, size N5x3
    extern double w5[7];      ///< Weights, size N5

    /// Order-17 (61-point) Dunavant rule on a triangle.
    extern int N17;
    extern double x17[61][3]; ///< Barycentric nodes, size N17x3
    extern double w17[61];    ///< Weights, size N17
} // namespace dunavant

/// \ingroup quadrature
namespace gausslegendre {
  /// Order-1 Gauss–Legendre rule on [0,1].
  extern int N1;
  extern double x1[1][2];  ///< @f$ (\xi, 1-\xi) @f$ nodes, size N1x2
  extern double w1[1];     ///< Weights, size N1

  /// Order-3 Gauss–Legendre rule on [0,1].
  extern int N3;
  extern double x3[3][2];  ///< @f$ (\xi, 1-\xi) @f$ nodes, size N3x2
  extern double w3[3];     ///< Weights, size N3

  /// Order-4 Gauss–Legendre rule on [0,1].
  extern int N4;
  extern double x4[4][2];  ///< @f$ (\xi, 1-\xi) @f$ nodes, size N4x2
  extern double w4[4];     ///< Weights, size N4

  /// Order-8 Gauss–Legendre rule on [0,1].
  extern int N8;
  extern double x8[8][2];  ///< @f$ (\xi, 1-\xi) @f$ nodes, size N8x2
  extern double w8[8];     ///< Weights, size N8

  /// Order-17 Gauss–Legendre rule on [0,1].
  extern int N17;
  extern double x17[68][2];///< @f$ (\xi, 1-\xi) @f$ nodes, size N17x2
  extern double w17[68];   ///< Weights, size N17
} // namespace gausslegendre
#endif
