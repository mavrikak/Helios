/**
 * @file Bsp3D.h
 * @brief 3D Binary Space Partitioning (BSP) utilities for triangle meshes.
 *
 * Declares a lightweight BSP tree used to spatially partition triangle meshes
 * and perform robust point‑location queries (inside/outside) with respect to
 * the closed surface defined by the mesh.
 */

#ifndef BSP_H
#define BSP_H

#include <cstdio>
#include <vector>

#include "Triangle.h"
#include "globals.h"

/**
 * @brief Polygon wrapper carried through BSP splitting.
 *
 * Each @c bspPolygon is initialized from a mesh triangle but stores its own
 * working list of nodes (rvec vertices) so that it can be split when the
 * polygon straddles a splitter plane. The field @c side indicates which face
 * orientation (0 or 1) defines the "outside" half‑space.
 */
struct bspPolygon {
  Triangle *origin;         //!< Source triangle.
  std::vector<rvec> nodes;  //!< Current polygon vertices.
  int side;                 //!< Orientation side (0 or 1).
};

/// \ingroup mesh
/**
 * @brief Minimal 3D BSP tree for triangle meshes.
 *
 * Usage:
 * - Add mesh elements as polygons via AddElement().
 * - Pick a splitter via SetSplitter() (last element in the list).
 * - Call Branch() to recursively partition elements into inside/outside
 * children until all leaves are marked as buds.
 * - Use TestPoint() to classify query points (inside or outside).
 */
class Bsp3D {
 private:
  bspPolygon *splitter;               //!< Polygon defining the splitter plane at this node.
  bool bud;                           //!< True if this node is a leaf (bud).
  Bsp3D *bspInside;                   //!< Child covering inside half‑space.
  Bsp3D *bspOutside;                  //!< Child covering outside half‑space.
  std::vector<bspPolygon *> elements; //!< Polygons queued for partitioning.

 public:
  /** @name Lifecycle */
  ///@{
  Bsp3D();  //!< Construct an empty node (non-bud by default).
  ~Bsp3D(); //!< Recursively deletes children (but not polygon origins).
  ///@}

  /** @name Build operations */
  ///@{
  /**
   * @brief Select the last pushed element as splitter at this node.
   *
   * Removes the chosen polygon from @c elements and stores its pointer in
   * @c splitter. The caller must ensure that @c elements is non-empty.
   * @return Always @c true.
   */
  bool SetSplitter();

  /**
   * @brief Add a triangle as a new polygon element.
   * @param inElement Triangle pointer.
   * @param inSide Face side (0 or 1) used to set the plane orientation.
   * @return Always @c true.
   */
  bool AddElement(Triangle *inElement, int inSide);

  /**
   * @brief Add an already constructed polygon to this node.
   * @param poly Polygon allocated by the caller.
   * @return Always @c true.
   */
  bool AddElement(bspPolygon *poly);

  /**
   * @brief Recursively branch on the current splitter and distribute elements.
   *
   * Creates @c bspInside and @c bspOutside children, assigns/splits every
   * pending polygon into the appropriate child list, then recurses. A child
   * receiving no elements becomes a bud.
   * @return Always @c true.
   */
  bool Branch();
  ///@}


  /** @name Node state */
  ///@{
  /**
   * @brief Mark this node as a leaf (bud).
   * @return Always @c true.
   */
  bool SetBud();

  /**
   * @brief Query whether this node is a bud (leaf).
   * @return @c true if leaf; @c false otherwise.
   */
  bool IsBud();
  ///@}


  /** @name Geometry queries */
  ///@{
  /**
   * @brief Classify a point relative to the splitter plane.
   * @param test Query position.
   * @return +1 if outside (positive distance), -1 if inside (negative
   * distance), 0 if coplanar to the splitter within a small tolerance.
   */
  int LocTest(rvec test);

  /**
   * @brief Signed distance from a point to the splitter plane.
   *
   * The plane is anchored at the splitter triangle barycenter and oriented by
   * its forward/backward normal according to @c splitter->side.
   * @param test Query position.
   * @return Signed distance (positive = outside half-space).
   */
  double LocDist(rvec test);

  /**
   * @brief Clip a polygon against the splitter, dispatching fragments to children.
   *
   * Evaluates all vertices against the splitter plane and, if the polygon
   * straddles the plane, inserts intersection points and produces two new
   * polygons: one for the inside child and one for the outside child. If the
   * polygon lies fully on one side, it is forwarded intact.
   *
   * @param poly Polygon to classify/split (input and, if fully on one side,
   * forwarded to a child).
   * @return @c true if @p poly was split into two newly allocated polygons and
   * the caller may delete the original; @c false if @p poly remains
   * owned/needed by a child.
   */
  bool SplitPoly(bspPolygon *poly);

  /**
   * @brief Recursive point-location query through the BSP tree.
   * @param point Query position.
   * @return +1 if the point ends in an outside bud; -1 if it ends in an inside
   * bud. Points coplanar to a splitter follow the @c inside branch by
   * convention.
   */
  int TestPoint(rvec point);
  ///@}
};

#endif
