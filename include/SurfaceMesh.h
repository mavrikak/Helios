/**
 * @file SurfaceMesh.h
 * @brief Declaration of the SurfaceMesh class and mesh edge struct.
 *
 * Holds nodes, triangles, and unique edges of a surface mesh; supports
 * reading/writing meshes (custom .mesh, COMSOL .mphtxt, Gmsh .msh),
 * periodic translations/symmetries, and a few utilities for geometry
 * queries used by RWG-based SIE solvers.
 */

// Only include this once during compiling
#ifndef SURFACEMESH_H
#define SURFACEMESH_H

#include <set>
#include <utility>
#include <vector>
#include "JobParser.h"
#include "Triangle.h"
#include "globals.h"
#include "mathFunctions.h"

/**
 * @struct edge
 * @brief Mesh edge record (unique and unoriented for membership; oriented for vectors).
 */
struct edge {
  rvec *node1;        ///< Pointer to 1st node within the mesh node array
  rvec *node2;        ///< Pointer to 2nd node within the mesh node array
  double edgeLength;  ///< cached Euclidean length |node2 - node1|
};

/// \ingroup mesh
/**
 * @class SurfaceMesh
 * @brief Lightweight surface mesh container with file I/O and BC helpers.
 *
 * Responsibilities:
 *  - Store nodes, triangles, and unique edges (each geometric edge appears once).
 *  - Read meshes from several formats (.mesh, .mphtxt, .msh) and write .mesh.
 *  - Provide edge queries (index/length/vector) and triangle lists per domain.
 *  - Handle periodic boundary conditions via stored symmetry operations.
 */
class SurfaceMesh {
 private:
  JobParser jobParser;                  ///< Internal parser for conversion/aux I/O
  std::vector<rvec> nodes;              ///< Mesh nodes (owned)
  std::vector<Triangle> triangles;      ///< Mesh triangles (owned)
  std::vector<edge> edges;              ///< Unique edges discovered from triangles
  std::set<int> domainIndeces;          ///< All domain IDs appearing on any triangle
  std::vector<rvec*> supNodes;          ///< Supplementary nodes for boundary conditions
  std::vector<Triangle*> supTriangles;  ///< Supplementary triangles for boundary conditions

  /**
   * @brief Insert edge (n1,n2) if not already present (directly or via BC image).
   * @param n1 First endpoint pointer.
   * @param n2 Second endpoint pointer.
   * @return 1 if a new edge is inserted; 0 if it already existed.
   */
  int AddEdgeTested(rvec* n1, rvec* n2);
  
  /**
   * @brief Return index of a node pointer within the owned node array.
   * @param in Pointer to a node inside @p nodes.
   * @return Zero-based index, or -1 if not found.
   */
  int GetNodeIndex(rvec* in);
  
  /// Boundary conditions (symmetry operations) used to detect periodic edges.
  std::vector<SymmetryOperation*> bcs;

 public:
  /** @brief Default constructor. */
  SurfaceMesh();
  
  /** @brief Destructor: frees BC ops and supplementary nodes/triangles. */
  virtual ~SurfaceMesh();
  
  /**
   * @brief Load a mesh from file by extension.
   *
   * Dispatches to:
   *  - LoadFromMphtxtFile() for ".mphtxt"
   *  - LoadFromMshFile()    for ".msh"
   *  - LoadFromMeshFile()   otherwise (custom .mesh)
   *
   * @param fileName Path to input mesh file.
   * @return 0 on success; non-zero on error.
   */
  int LoadFromFile(std::string fileName);
  
  /**
   * @brief Read mesh from the custom ".mesh" text format.
   * @param fileName Path to .mesh file.
   * @return 0 on success; non-zero on error.
   */
  int LoadFromMeshFile(std::string fileName);
  
  /**
   * @brief Read mesh from a COMSOL ".mphtxt" file.
   * @param fileName Path to .mphtxt file.
   * @return 0 on success; non-zero on error.
   */
  int LoadFromMphtxtFile(std::string fileName);
  
  /**
   * @brief Read mesh from a Gmsh ".msh" file (v2 text).
   * @param fileName Path to .msh file.
   * @return 0 on success; non-zero on error.
   */
  int LoadFromMshFile(std::string fileName);
  
  /**
   * @brief Write current mesh to the custom ".mesh" format.
   * @param fileName Output base name (extension not added automatically).
   * @return 0 on success; non-zero on error.
   */
  int SaveToMeshFile(std::string fileName);
  
  /**
   * @brief Highest domain index present in the mesh.
   * @return Max element of @p domainIndeces.
   */
  int MaxDomainIndex();
  
  /**
   * @brief Find the index of an existing edge between two nodes (direct match).
   * @param n1 First endpoint pointer.
   * @param n2 Second endpoint pointer.
   * @return Edge index, or -1 if not found.
   */
  int FindEdgeFB(rvec* n1, rvec* n2);
  
  /**
   * @brief Find the index of an edge mapped between n1,n2 via a stored symmetry op.
   * @param n1 First endpoint pointer.
   * @param n2 Second endpoint pointer.
   * @return {edgeIndex, symmetryOp*} if found; {-1, nullptr} otherwise.
   */
  std::pair<int, SymmetryOperation*> FindEdgeBC(rvec* n1, rvec* n2);
  
  /** 
   * @brief Return length of an edge by index. 
   * @param edgeIndex Edge index.
   * @return Edge length.
   */
  double EdgeLength(int edgeIndex);
  
  /** 
   * @brief Return oriented vector (node2 - node1) for edge index. 
   * @param edgeIndex Edge index.
   * @return Edge vector.
   */
  rvec EdgeVector(int edgeIndex);
  
  /** 
   * @brief Return a copy of the edge record (node pointers + length). 
   * @param edgeIndex Edge index.
   * @return edge struct.
   */
  edge Edge(int edgeIndex);
  
  /**
   * @brief Triangles bordering a given domain (FRONT or BACK).
   * @param dIndex Domain ID.
   * @return List of pointers to triangles touching @p dIndex.
   */
  std::vector<Triangle*> BorderingTriangles(int dIndex);
  
  /** 
   * @brief Total number of stored unique edges. 
   * @return Number of edges.
   */
  int EdgeCount();
  
  /**
   * @brief Register a new periodic translation (adds a symmetry operation).
   * @param translationVector Translation vector.
   * @return 0 on success.
   */
  int NewTranslation(rvec);
  
  /**
   * @brief Add a supplementary node (for boundary-condition images).
   * @param node Node position to store.
   * @return 0 on success.
   */
  int AddNode(rvec);
  
  /**
   * @brief Add a supplementary triangle (for boundary-condition images).
   * @param n1,n2,n3 Node pointers.
   * @param frontIndex,backIndex Bordering domain IDs.
   * @return 0 on success.
   */
  int AddTriangle(rvec* n1, rvec* n2, rvec* n3, int frontIndex, int backIndex);
  
  /** 
   * @brief Pointer to the last supplementary node added. 
   * @return Pointer to last node.
   */
  rvec* LastNode();
  
  /** 
   * @brief Pointer to the last supplementary triangle added. 
   * @return Pointer to last triangle.
   */
  Triangle* LastTriangle();
  
  /**
   * @brief Convert/modify a mesh according to a conversion job file.
   *
   * Recognized blocks:
   *  - @f$<mesh>@f$file@f$</mesh>@f$     : load a mesh
   *  - @f$<index>@f$n1 n2@f$</index>@f$  : SwitchIndex(n1,n2)
   *  - @f$<plane>...</plane>@f$          : RemovePlane(normal, point)
   *  - @f$<out>@f$base@f$</out>@f$       : SaveToMeshFile(base + ".mesh")
   *
   * @param inputMeshFile Path to conversion job file.
   * @return 0 on success.
   */
  int ConvertMeshFile(std::string inputMeshFile);
  
  /**
   * @brief Remove all triangles on or below a plane (strict > 1e-10 kept).
   * @param normal Plane normal.
   * @param point  A point on the plane.
   * @return 0 on success.
   */
  int RemovePlane(rvec normal, rvec point);
  
  /**
   * @brief Replace all occurrences of domain ID n1 with n2 in triangles.
   * @param n1 Old domain ID.
   * @param n2 New domain ID.
   * @return 0 on success.
   */
  int SwitchIndex(int n1, int n2);
  
  /**
   * @brief Minimum and maximum z among all nodes.
   * @return {zmin, zmax}; throws if mesh is empty.
   */
  std::pair<double, double> zMinMax();
};

#endif
