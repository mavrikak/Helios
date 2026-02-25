/**
 * @file RWGFun.h
 * @brief Header for the RWGFun class. The RWGFun class represents a Rao–Wilton–Glisson (RWG) basis function
 * defined on a single triangular facet (half RWG). It stores the owner triangle, the "free" vertex, 
 * the associated edge in the mesh, and whether the function is inverted w.r.t. the facet's FRONT/BACK side.
 */

// Only include this once during compiling
#ifndef RWGFUN_H
#define RWGFUN_H

#include "SurfaceMesh.h"
#include "Triangle.h"
#include "globals.h"

/// \ingroup mesh
/**
 * @class RWGFun
 * @brief RWG basis function tied to a mesh triangle.
 *
 * Given a triangle and one of its vertices designated as the "free" vertex,
 * the RWG function is (up to scale) linear along the triangle and proportional
 * to (r − r_free) (or the negative thereof) depending on orientation.
 *
 * This class offers:
 *  - Construction from mesh, triangle, free vertex and side (FRONT/BACK).
 *  - Pointwise evaluation of the vector shape function on the triangle.
 *  - Accessors for triangle/edge bookkeeping and sign/orientation.
 */
class RWGFun {
 private:
  Triangle* triangle; ///< Owning triangle for this half-RWG
  rvec* freeVertex;   ///< Pointer to the free vertex of the RWG on this triangle
  int edgeIndex;      ///< Index of the opposite edge in the mesh (across from freeVertex)
  bool inverted;      ///< True if function points towards free vertex (depends on side/orientation)
  SurfaceMesh* mesh;  ///< Mesh that contains the triangle and its edge
  int side;           ///< Triangle side used at construction (FRONT or BACK)

 public:
  /**
   * @brief Construct an RWG basis function tied to one triangle.
   *
   * Determines the opposite edge from @p inFreeVertex, resolves boundary
   * conditions (periodic/symmetry) via SurfaceMesh helpers if needed,
   * and sets the orientation (inverted) based on @p inSide and geometry.
   *
   * @param inMesh       Owning mesh.
   * @param inTri        Pointer to the triangle supporting this half-RWG.
   * @param inFreeVertex Pointer to the free vertex of the RWG on @p inTri.
   * @param inSide       FRONT or BACK (affects sign/orientation).
   */
  RWGFun(SurfaceMesh *inMesh, Triangle *inTri, rvec *inFreeVertex, int inSide);

  /**
   * @brief Evaluate the RWG vector shape function at a point on the triangle.
   *
   * Returns @f$\pm@f$(r − r_free) * L_edge / (2 * A_triangle), where the sign depends
   * on orientation (inverted flag). The caller must ensure @p r lies on the
   * triangle; no membership check is performed here.
   *
   * @param r Point of evaluation (in the triangle's plane).
   * @return RWG vector value at @p r.
   */
  rvec Evaluate(rvec r);
  
  /**
   *  @brief Get the pointer to the owner Triangle.
   *  @return Pointer to the owner Triangle.
   */
  Triangle *TrianglePtr();
  
  /** 
   * @brief Get the free vertex position (by value).
   * @return Free vertex position.
   */
  rvec FreeVertex();
  
  /** 
   * @brief Get the pointer to the free vertex.
   * @return Pointer to the free vertex.
   */
  rvec *FreeVertexPtr();
  
  /** 
   * @brief True if the RWG is inverted w.r.t. the triangle's FRONT side.
   * @return True if inverted.
   */
  bool Inverted();
  
  /**
   * @brief Return @f$\pm@f$L_n / A_n (missing factor 1/2 by convention).
   * @return RWG pre-factor.
   *
   * This equals the surface divergence of the RWG function on the triangle.
   * Sign depends on @ref Inverted().
   */
  double RWGPreFactor();
  
  /** 
   * @brief Index of the opposite edge in the mesh.
   * @return Index of the opposite edge in the mesh.
   */
  int EdgeIndex();
  
  /** 
   * @brief Return the opposite edge structure from the mesh.
   * @return Opposite edge structure from the mesh.
   */
  edge Edge();

  /** 
   * @brief Stored side (FRONT or BACK) used at construction.
   * @return Stored side.
   */
  int Side();
};

#endif

