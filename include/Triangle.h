/**
 * @file Triangle.h
 * @brief Declaration of the Triangle class (surface element with cached geometry).
 *
 * The Triangle class stores pointers to its three corner nodes and the two
 * bordering domain indices (front/back). It caches geometric quantities
 * (reference normal, area, perimeter, center of gravity) and provides
 * analytic self-integrals.
 */

#ifndef TRIANGLE_H
#define TRIANGLE_H

#include <vector>
#include "globals.h"

/// \ingroup mesh
/**
 * @class Triangle
 * @brief Triangular surface element with cached geometry and utilities.
 *
 * Stores node pointers and bordering domains, computes and caches:
 *  - reference unit normal (pointing to FRONT),
 *  - surface area and perimeter,
 *  - centroid (center of gravity),
 *  - analytic integrals I1, I2s, I2d.
 *
 * Notes:
 *  - Node accessors return by value (Node) or pointer (NodePtr).
 *  - FBNormalVector() flips the reference normal for BACK.
 */
class Triangle {
 private:
  // --- data ---
  std::vector<rvec *> nodes;      ///< Pointers to the three corner nodes
  std::vector<int> domainIndeces; ///< [FRONT, BACK] domain indices
  rvec normalVector;              ///< Unit normal pointing to FRONT
  double surfaceArea;             ///< Triangle area
  double perimeter;               ///< Triangle perimeter (sum of edge lengths)
  double I1;                      ///< Analytic self integral (R^2 over triangle)
  double I2s[3];                  ///< Analytic self integrals (single-vertex cases)
  double I2d[3];                  ///< Analytic self integrals (double-vertex cases)
  rvec center;                    ///< Centroid
  
 public:
  // ---------------- Constructors & basic setup ----------------
  /**
   * @brief Construct from 3 node pointers and bordering domains.
   * @param n1,n2,n3  Node pointers.
   * @param dIndex1   Front domain index.
   * @param dIndex2   Back  domain index.
   */
  Triangle(rvec *n1, rvec *n2, rvec *n3, int dIndex1, int dIndex2);
  
  /**
   * @brief Construct from 3 node pointers; domains can be set later.
   * @param n1,n2,n3  Node pointers.
   */
  Triangle(rvec *n1, rvec *n2, rvec *n3);
  
  /// @brief Default constructor (uninitialized geometry).
  Triangle(){};
  
  /**
   * @brief Replace the three node pointers and recompute cached geometry.
   * @param n1,n2,n3  Node pointers.
   * @return 0 on success.
   */
  int AssignNodes(rvec *n1, rvec *n2, rvec *n3);

  /**
   * @brief Set front/back domain indices; recomputes cached geometry.
   * @param dIndex1 Front domain index.
   * @param dIndex2 Back  domain index.
   * @return 0 on success.
   */
  int SetDomainIndeces(int dIndex1, int dIndex2);
  
  /**
   * @brief Compute/cached geometry (normal, area, perimeter, center, I1/I2s/I2d).
   * @return 0 on success.
   */
  int CalculateGeometryData();
  
  /**
   * @brief Invert orientation (swap first two nodes).
   * @return 0 on success.
   */
  int Invert();
  
  // ---------------- Accessors ----------------
  /** 
   * @brief Node position (by value) for local index 0..2.
   * @param index Local index.
   * @return Node position.
   */
  rvec Node(int index);
  
  /** 
   * @brief Node pointer for local index 0..2.
   * @param index Local index.
   * @return Node pointer.
   */
  rvec *NodePtr(int index);
  
  /**
   * @brief Whether a domain id matches FRONT/BACK; returns FRONT, BACK, or −1.
   * @param dIndex Domain index.
   * @return FRONT, BACK, or −1.
   */
  int DomainPosition(int dIndex);
  
  /** 
   * @brief Reference unit normal (points to FRONT).
   * @return Reference normal.
   */
  rvec RefNormalVector();
  
  /**
   * @brief FRONT/BACK normal: FRONT -> +normal, BACK -> −normal.
   * @param side FRONT or BACK.
   * @return Normal.
   */
  rvec FBNormalVector(int side);
  
  /** 
   * @brief Surface area of the triangle.
   * @return Area.
   */
  double Area();

  /** 
   * @brief Perimeter (sum of edge lengths).
   * @return Perimeter.
   */
  double Perimeter();

  /** 
   * @brief Centroid of the triangle.
   * @return Centroid.
   */
  rvec Center();
  
  /** 
   * @brief Bordering domain indices as {front, back}.
   * @return {front, back}.
   */
  std::pair<int, int> BorderingIndeces();

  /**
   * @brief Double integral of R^2 over the triangle (analytic, self term).
   * @return Double integral.
   */
  double IntRR();

  /**
   * @brief Replace a bordering domain index n1 with n2 (if present).
   * @param n1,n2 Bordering domain indices.
   * @return 0 on success.
   */
  int SwitchIndex(int n1, int n2);
  
  // ---------------- Adjacency / common entities ----------------
  /**
   * @brief Count common nodes with triangle T (0..3).
   * @param T Triangle.
   * @return Number of shared nodes.
   */
  int FindAdjacency(Triangle *T);
  
  /**
   * @brief Find common edge local indices with T.
   * @param T Triangle.
   * @return {i+1, -(j+1)} if edges match (sign conveys orientation), or {0,0}.
   */
  std::pair<int, int> FindCommonEdgeIndices(Triangle *T);
  
  /**
   * @brief Find common vertex local indices with T.
   * @param T Triangle.
   * @return {i+1, j+1} if a shared node exists, or {0,0}.
   */
  std::pair<int, int> FindCommonVertexIndices(Triangle *T);
  
  /**
   * @brief Count common nodes with T after translating T by vector t.
   * @param T Triangle.
   * @param t Translation vector.
   * @return Number of shared nodes.
   */
  int FindAdjacency(Triangle *T, rvec t);

  // ---------------- Point-in-triangle test ----------------
  /**
   * @brief Test whether a point lies on (or very near) the triangle.
   *
   * The point is projected to the triangle plane and tested using
   * barycentric/simplex coordinates with tolerance @p tol.
   *
   * @param p   Observation point.
   * @param tol Tolerance for coplanarity and barycentric tests.
   * @return true if p is on/near the facet; false otherwise.
   */
  bool containsPoint(const rvec &p, double tol = 1e-6);

  /** 
   * @brief Small helper to zero tiny values (for robust barycentric math).
   * @param x Value.
   * @param eps Tolerance.
   * @return Zeroed value.
   */
  double clean(double x, double eps = 1e-10) {
    return std::abs(x) < eps ? 0.0 : x;
  }

  // ---------------- Analytic self-integrals ----------------
  /** 
   * @brief Return cached I1.
   * @return I1.
   */
  double getI1();

  /** 
   * @brief Return cached I2 (single-vertex) for node pointer r.
   * @param r Node pointer.
   * @return I2.
   */
  double getI2(rvec *r);

  /** 
   * @brief Return cached I2 (double-vertex) for node pair (r1,r2).
   * @param r1,r2 Node pointers.
   * @return I2.
   */
  double getI2(rvec *r1, rvec *r2);

  /**
   * @brief Return I2 for a translated configuration (periodic images).
   * @param r1,r2 Endpoints in original cell.
   * @param t Translation vector for r1.
   * @return I2.
   */
  double getI2Translated(rvec r1, rvec r2, rvec t);
};

#endif

