/**
 * @file RWGFun.cpp
 * @brief Implementation for the RWGFun class.
 */

#include "RWGFun.h"

/**
 * @details Locates the edge opposite to @p inFreeVertex, handles boundary-condition
 * edges (periodic/symmetry) by mapping and augmenting the mesh as needed,
 * and determines the sign (inverted) based on the requested @p inSide and
 * the geometric orientation between edge/normal/free-vertex.
 */
RWGFun::RWGFun(SurfaceMesh *inMesh, Triangle *inTri, rvec *inFreeVertex,
               int inSide)
    : triangle(inTri), freeVertex(inFreeVertex), mesh(inMesh), side(inSide) {
  rvec *n1, *n2;

  // Determine edge vertices
  if (triangle->NodePtr(0) == freeVertex) {
    n1 = triangle->NodePtr(1);
    n2 = triangle->NodePtr(2);
  } else if (triangle->NodePtr(1) == freeVertex) {
    n1 = triangle->NodePtr(0);
    n2 = triangle->NodePtr(2);
  } else {
    n1 = triangle->NodePtr(0);
    n2 = triangle->NodePtr(1);
  }

  // Determine edge index
  edgeIndex = mesh->FindEdgeFB(n1, n2);

  // If any boundary conditions on the mesh
  if (edgeIndex == -1) {
    std::pair<int, SymmetryOperation *> result(mesh->FindEdgeBC(n1, n2));
    SymmetryOperation *s(result.second);
    
    // Map the free vertex across the symmetry and append as a new node
    mesh->AddNode(s->Apply(*freeVertex));
    freeVertex = mesh->LastNode();
    edgeIndex = result.first;

    // Ensure (n1, n2) match the edge orientation just resolved
    std::pair<int, int> d(triangle->BorderingIndeces());
    double eps(1.e-2);
    rvec dif(s->Apply(*n1) - *mesh->Edge(edgeIndex).node1);
    if (dot(dif, dif) < eps) {
      n1 = mesh->Edge(edgeIndex).node1;
      n2 = mesh->Edge(edgeIndex).node2;
    } else {
      n1 = mesh->Edge(edgeIndex).node2;
      n2 = mesh->Edge(edgeIndex).node1;
    }

    // Create the mirrored triangle with consistent orientation
    rvec normalVector(cross((rvec)(*n2 - *n1), (rvec)(*freeVertex - *n1)));
    if (dot(normalVector, triangle->FBNormalVector(FRONT)) > 0) {
      mesh->AddTriangle(n1, n2, freeVertex, d.first, d.second);
    } else {
      mesh->AddTriangle(n2, n1, freeVertex, d.first, d.second);
    }
    triangle = mesh->LastTriangle();
  }

  // Determine sign of RWG function
  if (inSide == FRONT)
    inverted = false;
  else
    inverted = true;

  // Adjust orientation using geometric consistency
  if (dot(cross(mesh->EdgeVector(edgeIndex), triangle->RefNormalVector()),
          (rvec)(*n1 - *freeVertex)) < 0)
    inverted = !inverted;
}

/**
 * @details Evaluate the RWG vector function at point @p r on the triangle.
 *
 * Uses the standard half-RWG formula:
 *  f(r) = @f$\pm@f$ (r − r_free) * (L_edge) / (2 A_tri),
 * where the sign is controlled by @ref Inverted().
 */
rvec RWGFun::Evaluate(rvec r) {
  if (inverted)
    return (*freeVertex - r) * mesh->EdgeLength(edgeIndex) /
           (2. * triangle->Area());
  else
    return (r - *freeVertex) * mesh->EdgeLength(edgeIndex) /
           (2. * triangle->Area());
}

// Owner triangle pointer.
Triangle *RWGFun::TrianglePtr() { return triangle; }

// Free vertex (by value).
rvec RWGFun::FreeVertex() { return *freeVertex; }

// Free vertex pointer.
rvec *RWGFun::FreeVertexPtr() { return freeVertex; }

// True if oriented towards the free vertex.
bool RWGFun::Inverted() { return inverted; }

/**
 * @details Return @f$\pm@f$L_n / A_n (missing factor 1/2), 
 * equal to surface divergence. Sign reflects @ref Inverted().
 */
double RWGFun::RWGPreFactor() {
  if (inverted)
    return -mesh->EdgeLength(edgeIndex) / triangle->Area();
  else
    return mesh->EdgeLength(edgeIndex) / triangle->Area();
}

// Mesh edge index opposite the free vertex.
int RWGFun::EdgeIndex() { return edgeIndex; }

// Mesh edge structure opposite the free vertex.
edge RWGFun::Edge() { return mesh->Edge(edgeIndex); }

// Stored side passed at construction (FRONT/BACK).
int RWGFun::Side() { return side; }
