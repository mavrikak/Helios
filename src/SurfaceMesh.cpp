/**
 * @file SurfaceMesh.cpp
 * @brief Implementation of SurfaceMesh mesh I/O and utilities.
 */

#include "SurfaceMesh.h"
#include <algorithm>
#include <fstream>
#include <iostream>

#include "iofunctions.h"

// Default constructor.
SurfaceMesh::SurfaceMesh() {}

/**
 * @details Destructor: free symmetry ops and supplementary mesh entities.
 *
 * Frees:
 *  - all SymmetryOperation* stored in @p bcs
 *  - supplementary nodes added via AddNode()
 *  - supplementary triangles added via AddTriangle()
 */
SurfaceMesh::~SurfaceMesh() {
  for (unsigned int i(0); i < bcs.size(); ++i) {
    delete bcs[i];
  }
  for (unsigned int i(0); i < supNodes.size(); ++i) {
    delete supNodes[i];
  }
  for (unsigned int i(0); i < supTriangles.size(); ++i) {
    delete supTriangles[i];
  }
}

// Load mesh by extension: .mphtxt -> COMSOL, .msh -> Gmsh, else .mesh.
int SurfaceMesh::LoadFromFile(std::string fileName) {
  supNodes.clear();
  supTriangles.clear();
  if (fileName.length() > 6)
    if (fileName.substr(fileName.length() - 7).compare(".mphtxt") == 0)
      return LoadFromMphtxtFile(fileName);
  if (fileName.length() > 3)
    if (fileName.substr(fileName.length() - 4).compare(".msh") == 0) {
      return LoadFromMshFile(fileName);
    }
  return LoadFromMeshFile(fileName);
}

/**
 * @details Read custom ".mesh" format: nodes, then triangles (+ domains), detect edges.
 *
 * Expected layout:
 *   @f$<Nnodes>@f$
 *   x y z   (repeated Nnodes times)
 *   @f$<Ntris>@f$
 *   n1 n2 n3 d1 d2   (1-based node indices, two bordering domains)
 */
int SurfaceMesh::LoadFromMeshFile(std::string fileName) {
  // Open mesh file
  std::ifstream fin(fileName.c_str());
  if (!fin.is_open()) {
    std::cout << "Error opening mesh file: " << fileName << std::endl;
    return 1;
  }
  int cnt;
  double nx, ny, nz;

  // Read number of nodes
  int err = ReadCommented<int>(&fin, &cnt);
  if (err) {
    std::cout << "Error reading mesh file: " << fileName << std::endl;
    return 1;
  }
  // Read nodes
  nodes.clear();
  for (int i = 0; i < cnt; i++) {
    err = ReadCommented<double>(&fin, &nx) + ReadCommented<double>(&fin, &ny) +
          ReadCommented<double>(&fin, &nz);
    if (err) {
      std::cout << "Error reading mesh file: " << fileName << std::endl;
      return 1;
    }
    nodes.push_back(rvec(nx, ny, nz));
  }
  std::cout << "Read " << cnt << " nodes." << std::endl;

  int n1, n2, n3, d1, d2;

  // Read number of triangles
  err = ReadCommented<int>(&fin, &cnt);
  if (err) {
    std::cout << "Error reading mesh file: " << fileName << std::endl;
    return 1;
  }
  // Read triangles and detect edges
  triangles.clear();
  domainIndeces.clear();
  for (int i = 0; i < cnt; i++) {
    err = ReadCommented<int>(&fin, &n1) + ReadCommented<int>(&fin, &n2) +
          ReadCommented<int>(&fin, &n3) + ReadCommented<int>(&fin, &d1) +
          ReadCommented<int>(&fin, &d2);
    if (err) {
      std::cout << "Error reading mesh file: " << fileName << std::endl;
      return 1;
    }
    triangles.push_back(
        Triangle(&(nodes[n1 - 1]), &(nodes[n2 - 1]), &(nodes[n3 - 1]), d1, d2));
    domainIndeces.insert(d1);
    domainIndeces.insert(d2);

    // Detect edges
    AddEdgeTested(&(nodes[n1 - 1]), &(nodes[n2 - 1]));
    AddEdgeTested(&(nodes[n2 - 1]), &(nodes[n3 - 1]));
    AddEdgeTested(&(nodes[n3 - 1]), &(nodes[n1 - 1]));
  }
  std::cout << "Read " << triangles.size() << " triangles and found "
            << edges.size() << " edges." << std::endl;
  fin.close();
  return 0;
}

/**
 * @details Read mesh from a COMSOL 4.2 ".mphtxt" file and border-domain block.
 *
 * Reads the node block, then the triangle connectivity (node indices),
 * then the bordering domain pairs (d1,d2). Edges are detected uniquely.
 */
int SurfaceMesh::LoadFromMphtxtFile(std::string fileName) {
  // Open mesh file
  std::ifstream fin(fileName.c_str());
  if (!fin.is_open()) {
    std::cout << "Error opening mesh file: " << fileName << std::endl;
    return 1;
  }

  int nBuf(0), err(0);
  std::string strBuf("");

  // Version numbers and no. of tags
  for (int i = 0; i < 3; i++) ReadCommented<int>(&fin, &nBuf);

  // Tags
  for (int i = 0; i < nBuf; i++) {
    ReadCommented<std::string>(&fin, &strBuf);
    getline(fin, strBuf);
  }

  // Number of types
  ReadCommented<int>(&fin, &nBuf);

  // Types
  for (int i = 0; i < nBuf; i++) {
    ReadCommented<std::string>(&fin, &strBuf);
    getline(fin, strBuf);
  }

  // Object header
  for (int i = 0; i < 3; i++) {
    ReadCommented<int>(&fin, &nBuf);
  }
  ReadCommented<int>(&fin, &nBuf);
  getline(fin, strBuf);

  // Version
  int version;
  ReadCommented<int>(&fin, &version);
  // Number of nodes and version
  int cnt(0);
  for (int i = 0; i < 2; i++) ReadCommented<int>(&fin, &cnt);
  // Lowest node index
  int lowestNode(0);
  ReadCommented<int>(&fin, &lowestNode);

  // Nodes
  nodes.clear();
  for (int i = 0; i < cnt; i++) {
    double nx, ny, nz;
    err = ReadCommented<double>(&fin, &nx) + ReadCommented<double>(&fin, &ny) +
          ReadCommented<double>(&fin, &nz);
    if (err) {
      std::cout << "Error reading mesh file: " << fileName << std::endl;
      return 1;
    }
    nodes.push_back(rvec(nx, ny, nz));
  }
  std::cout << "Read " << cnt << " nodes." << std::endl;

  // Jump to triangles
  err = 0;
  while (strBuf.compare("tri") != 0)
    err += ReadCommented<std::string>(&fin, &strBuf);

  // Number of triangles
  for (int i = 0; i < 2; i++) err += ReadCommented<int>(&fin, &cnt);
  if (err) {
    std::cout << "Error reading mesh file: " << fileName << std::endl;
    return 1;
  }

  // Read triangles
  triangles.clear();
  for (int i = 0; i < cnt; i++) {
    int n1, n2, n3;
    err = ReadCommented<int>(&fin, &n1) + ReadCommented<int>(&fin, &n2) +
          ReadCommented<int>(&fin, &n3);
    if (err) {
      std::cout << "Error reading mesh file: " << fileName << std::endl;
      return 1;
    }
    triangles.push_back(Triangle(&(nodes[n1 - lowestNode]),
                                 &(nodes[n2 - lowestNode]),
                                 &(nodes[n3 - lowestNode])));
  }

  // Read and discard parameters
  int nTmp;
  err = ReadCommented(&fin, &nTmp);
  if (version == 2)  // Accomodation to version 4.2 of Comsol
    nTmp *= 3;
  err += ReadCommented(&fin, &cnt);
  if (err) {
    std::cout << "Error reading mesh file: " << fileName << std::endl;
    return 1;
  }

  double fTmp;
  for (int i = 0; i < cnt * nTmp; i++) {
    err += ReadCommented<double>(&fin, &fTmp);
  }

  // Read and discard domains
  err = ReadCommented(&fin, &cnt);
  for (int i = 0; i < cnt; i++) err += ReadCommented<int>(&fin, &nTmp);

  // Read and process bordering domains
  err = ReadCommented(&fin, &cnt);
  domainIndeces.clear();
  edges.clear();
  for (int i = 0; i < cnt; i++) {
    int d1, d2;
    err += ReadCommented(&fin, &d1);
    err += ReadCommented(&fin, &d2);
    domainIndeces.insert(d1);
    domainIndeces.insert(d2);
    if (d1 > d2) {
      triangles.at(i).Invert();
      triangles.at(i).SetDomainIndeces(d2, d1);
    } else
      triangles.at(i).SetDomainIndeces(d1, d2);

    // Detect edges
    AddEdgeTested(triangles.at(i).NodePtr(0), triangles.at(i).NodePtr(1));
    AddEdgeTested(triangles.at(i).NodePtr(1), triangles.at(i).NodePtr(2));
    AddEdgeTested(triangles.at(i).NodePtr(2), triangles.at(i).NodePtr(0));
  }
  std::cout << "Read " << triangles.size() << " triangles and found "
            << edges.size() << " edges." << std::endl;
  fin.close();

  return 0;
}

/**
 * @details Read mesh from a Gmsh ".msh" file and an auxiliary conversion job.
 *
 * Uses JobParser to optionally map element IDs to (dom1,dom2) pairs via a
 * @f$<boundary>@f$ block inside the currently open job file (external to the .msh).
 */
int SurfaceMesh::LoadFromMshFile(std::string fileName) {
  // Open mesh file
  std::ifstream fin(fileName.c_str());
  if (!fin.is_open()) {
    std::cout << "Error opening mesh file: " << fileName << std::endl;
    return 1;
  }

  double dBuf(0);
  std::string strBuf("");

  // Skipping MeshFormat data
  for (int i = 0; i < 4; i++) ReadCommented<double>(&fin, &dBuf);
  int cnt(0);
  int err(0);

  // Read number of nodes
  for (int i = 0; i < 3; i++) ReadCommented<int>(&fin, &cnt);

  if (err) {
    std::cout << "Error reading mesh file: " << fileName << std::endl;
    return 1;
  }

  // Nodes
  nodes.clear();
  for (int i = 0; i < cnt; i++) {
    int index;
    double nx, ny, nz;
    err = ReadCommented<int>(&fin, &index) + ReadCommented<double>(&fin, &nx) +
          ReadCommented<double>(&fin, &ny) + ReadCommented<double>(&fin, &nz);
    if (err) {
      std::cout << "Error reading mesh file: " << fileName << std::endl;
      return 1;
    }
    nodes.push_back(rvec(nx, ny, nz));
  }
  std::cout << "Read " << cnt << " nodes." << std::endl;

  // Read number of elements
  for (int i = 0; i < 3; i++) err = ReadCommented<int>(&fin, &cnt);
  if (err) {
    std::cout << "Error reading mesh file: " << fileName << std::endl;
    return 1;
  }

  // Prepares domains assignment to triangles
  jobParser.readTag();
  jobParser.readTag();
  std::vector<std::vector<int> > indexes;
  if (jobParser.getTag() == "boundary") {
    while (jobParser.getTag() == "boundary") {
      std::vector<int> buffer;
      for (unsigned int i(0); i < 3; ++i)
        buffer.push_back(jobParser.read<int>());
      indexes.push_back(buffer);
      jobParser.readTag();
      jobParser.readTag();
    }
  }

  // Triangles and domain indeces
  triangles.clear();
  domainIndeces.clear();
  for (int i = 0; i < cnt; i++) {
    int index;
    int type;
    int elementary;
    int n1, n2, n3;
    int id1;
    int id2;
    err = ReadCommented<int>(&fin, &index) + ReadCommented<int>(&fin, &type);
    if (err) {
      std::cout << "Error reading mesh file: " << fileName << std::endl;
      return 1;
    }
    if (type == 15) {
      for (int i = 0; i < 4; i++) err = ReadCommented<int>(&fin, &elementary);
    } else if (type == 1) {
      for (int i = 0; i < 5; i++) err = ReadCommented<int>(&fin, &elementary);
    } else if (type == 2) {
      err = ReadCommented<int>(&fin, &elementary) +
            ReadCommented<int>(&fin, &id1) + ReadCommented<int>(&fin, &id2) +
            ReadCommented<int>(&fin, &n1) + ReadCommented<int>(&fin, &n2) +
            ReadCommented<int>(&fin, &n3);
      int dom1(0);
      int dom2(1);

      for (unsigned int j = 0; j < indexes.size(); ++j) {
        if (indexes[j][0] == id2) {
          dom1 = indexes[j][1];
          dom2 = indexes[j][2];
        }
      }

      triangles.push_back(Triangle(&(nodes[n1 - 1]), &(nodes[n2 - 1]),
                                   &(nodes[n3 - 1]), dom1, dom2));
      domainIndeces.insert(dom1);
      domainIndeces.insert(dom2);
      // Detect edges
      AddEdgeTested(&(nodes[n1 - 1]), &(nodes[n2 - 1]));
      AddEdgeTested(&(nodes[n2 - 1]), &(nodes[n3 - 1]));
      AddEdgeTested(&(nodes[n3 - 1]), &(nodes[n1 - 1]));
    } else {
      std::cout << "Error: wrong element type in mesh " << std::endl;
      return 1;
    }
  }
  std::cout << "Read " << triangles.size() << " triangles and found "
            << edges.size() << " edges." << std::endl;
  fin.close();
  return 0;
}

// Return the highest domain number present.
int SurfaceMesh::MaxDomainIndex() { return *(domainIndeces.end()--); }

// Find an existing edge by matching endpoint pointers exactly.
int SurfaceMesh::FindEdgeFB(rvec* n1, rvec* n2) {
  std::vector<edge>::iterator edgeIter;
  int pos = 0;
  for (edgeIter = edges.begin(); edgeIter != edges.end(); edgeIter++) {
    if ((edgeIter->node1 == n1 && edgeIter->node2 == n2) ||
        (edgeIter->node1 == n2 && edgeIter->node2 == n1))
      return pos;
    pos++;
  }
  return -1;
}

// Find edge separated by a stored symmetry operation (periodic image).
std::pair<int, SymmetryOperation*> SurfaceMesh::FindEdgeBC(rvec* n1, rvec* n2) {
  std::vector<edge>::iterator edgeIter;
  double eps(1.e-2);
  int pos = 0;
  for (edgeIter = edges.begin(); edgeIter != edges.end(); edgeIter++) {
    for (unsigned int i(0); i < bcs.size(); ++i) {
      rvec dif11(bcs[i]->Apply(*n1 - *(edgeIter->node1)));
      rvec dif12(bcs[i]->Apply(*n1 - *(edgeIter->node2)));
      rvec dif21(bcs[i]->Apply(*n2 - *(edgeIter->node1)));
      rvec dif22(bcs[i]->Apply(*n2 - *(edgeIter->node2)));
      if (((dot(dif11, dif11) < eps) && (dot(dif22, dif22) < eps)) ||
          ((dot(dif21, dif21) < eps) && (dot(dif12, dif12) < eps))) {
        std::pair<int, SymmetryOperation*> out(pos, bcs[i]);
        return out;
      }
    }
    pos++;
  }
  std::pair<int, SymmetryOperation*> out(-1, 0);
  return out;
}

// Write current mesh to the custom ".mesh" text format.
int SurfaceMesh::SaveToMeshFile(std::string fileName) {
  // Open mesh file
  std::ofstream fout(fileName.c_str());
  if (!fout.is_open()) {
    std::cout << "Error opening mesh output file: " << fileName << std::endl;
    return 1;
  }
  fout << nodes.size() << std::endl;
  fout << std::setprecision(8) << std::fixed;
  for (std::vector<rvec>::iterator nodIter = nodes.begin();
       nodIter != nodes.end(); nodIter++) {
    fout << (*nodIter)(0) << "\t" << (*nodIter)(1) << "\t" << (*nodIter)(2)
         << std::endl;
  }
  fout << triangles.size() << std::endl;
  for (std::vector<Triangle>::iterator triIter = triangles.begin();
       triIter != triangles.end(); triIter++) {
    int n1(GetNodeIndex(triIter->NodePtr(0)) + 1);
    int n2(GetNodeIndex(triIter->NodePtr(1)) + 1);
    int n3(GetNodeIndex(triIter->NodePtr(2)) + 1);
    int d1(triIter->BorderingIndeces().first);
    int d2(triIter->BorderingIndeces().second);
    if (n1 == -1 || n2 == -1 || n3 == -1) return 1;
    fout << n1 << "\t" << n2 << "\t" << n3 << "\t" << d1 << "\t" << d2
         << std::endl;
  }
  fout.close();
  return 0;
}

// Locate the index of a node pointer in the node array, or -1.
int SurfaceMesh::GetNodeIndex(rvec* in) {
  for (unsigned int i = 0; i < nodes.size(); i++) {
    if (in == &(nodes.at(i))) return i;
  }
  return -1;
}

// Edge length by index (cached).
double SurfaceMesh::EdgeLength(int edgeIndex) {
  return edges[edgeIndex].edgeLength;
}

// Return a copy of the edge record (node pointers + length).
edge SurfaceMesh::Edge(int edgeIndex) { return edges[edgeIndex]; }

// Oriented edge vector = node2 − node1.
rvec SurfaceMesh::EdgeVector(int edgeIndex) {
  return *(edges[edgeIndex].node2) - *(edges[edgeIndex].node1);
}

// Return triangles bordering a given domain (FRONT or BACK).
std::vector<Triangle*> SurfaceMesh::BorderingTriangles(int dIndex) {
  std::vector<Triangle*> list;
  for (unsigned i = 0; i < triangles.size(); i++) {
    if (triangles[i].DomainPosition(dIndex) == FRONT ||
        triangles[i].DomainPosition(dIndex) == BACK)
      list.push_back(&(triangles[i]));
  }
  return list;
}

// Total number of unique edges stored.
int SurfaceMesh::EdgeCount() { return edges.size(); }

// Add an edge to list if not already present
int SurfaceMesh::AddEdgeTested(rvec* n1, rvec* n2) {
  edge edgeTemp;
  if ((FindEdgeFB(n1, n2) == -1) && (FindEdgeBC(n1, n2).first == -1)) {
    edgeTemp.node1 = n1;
    edgeTemp.node2 = n2;
    rvec vTemp = *n1 - *n2;
    edgeTemp.edgeLength = sqrt(dot(vTemp, vTemp));
    edges.push_back(edgeTemp);
    return 1;
  }
  return 0;
}

// Add a periodic translation to the list of symmetry operations.
int SurfaceMesh::NewTranslation(rvec translationVector) {
  bcs.push_back(new Translation(translationVector));
  return 0;
}

// Add a supplementary node (owned).
int SurfaceMesh::AddNode(rvec node) {
  supNodes.push_back(new rvec(node[0], node[1], node[2]));
  return 0;
}

// Pointer to last supplementary node added.
rvec* SurfaceMesh::LastNode() { return supNodes[supNodes.size() - 1]; }

// Add a supplementary triangle (owned).
int SurfaceMesh::AddTriangle(rvec* n1, rvec* n2, rvec* n3, int d1, int d2) {
  supTriangles.push_back(new Triangle(n1, n2, n3, d1, d2));
  return 0;
}

// Pointer to last supplementary triangle added.
Triangle* SurfaceMesh::LastTriangle() {
  return supTriangles[supTriangles.size() - 1];
}

/**
 * @details Convert/modify a mesh according to a conversion job file.
 *
 * Recognized directives:
 *  - @f$<mesh>@f$file@f$</mesh>@f$    : LoadFromFile(file)
 *  - @f$<index>@f$n1 n2@f$</index>@f$ : SwitchIndex(n1, n2)
 *  - @f$<plane>@f$...@f$</plane>@f$   : RemovePlane(normal, point)
 *  - @f$<out>@f$base@f$</out>@f$      : SaveToMeshFile(base + ".mesh")
 */
int SurfaceMesh::ConvertMeshFile(std::string inputMeshFile) {
  jobParser.open(inputMeshFile);
  jobParser.readTag();
  if (jobParser.getTag() == "conversion") {
    while (jobParser.lastRead() != "/conversion") {
      jobParser.readTag();
      if (jobParser.getTag() == "mesh") {
        LoadFromFile(jobParser.read<std::string>());
      }
      if (jobParser.getTag() == "index") {
        int n1(jobParser.read<int>());
        int n2(jobParser.read<int>());
        SwitchIndex(n1, n2);
        jobParser.readTag();
      }
      if (jobParser.getTag() == "plane") {
        rvec normal(0, 0, 0);
        rvec point(0, 0, 0);
        while (jobParser.lastRead() != "/plane") {
          jobParser.readTag();
          if (jobParser.getTag() == "normal") {
            double x(jobParser.read<double>());
            double y(jobParser.read<double>());
            double z(jobParser.read<double>());
            normal = rvec(x, y, z);
            jobParser.readTag();
          }
          if (jobParser.getTag() == "point") {
            double x(jobParser.read<double>());
            double y(jobParser.read<double>());
            double z(jobParser.read<double>());
            point = rvec(x, y, z);
            jobParser.readTag();
          }
        }
        RemovePlane(normal, point);
      }
      if (jobParser.getTag() == "out") {
        std::cout << "After conversion, mesh has " << triangles.size()
                  << " triangles." << std::endl;
        SaveToMeshFile(jobParser.read<std::string>() + ".mesh");
        jobParser.readTag();
      }
    }
  }
  jobParser.close();
  return 0;
}

// Remove triangles at/below a plane: keep those with (center - point)·normal > 1e-10.
int SurfaceMesh::RemovePlane(rvec normal, rvec point) {
  std::vector<Triangle> newTriangles;
  for (unsigned int i(0); i < triangles.size(); ++i) {
    double prod(dot(triangles[i].Center() - point, normal));
    if (prod > 1e-10) newTriangles.push_back(triangles[i]);
  }
  triangles = newTriangles;
  return 0;
}

// Replace all occurrences of domain ID n1 with n2 across triangles.
int SurfaceMesh::SwitchIndex(int n1, int n2) {
  for (unsigned int i(0); i < triangles.size(); ++i) {
    triangles[i].SwitchIndex(n1, n2);
  }
  std::cout << "Domain index " << n1 << " replaced by " << n2 << "."
            << std::endl;
  return 0;
}

// Compute min and max z among nodes; throws if mesh has no nodes.
std::pair<double, double> SurfaceMesh::zMinMax() {
  if (nodes.empty())
    throw std::runtime_error("zMinMax: nodes is empty");

  double zmin = nodes[0][2];
  double zmax = zmin;

  for (size_t i = 1; i < nodes.size(); ++i) {
    const double z = nodes[i][2];
    if (z < zmin) zmin = z;
    if (z > zmax) zmax = z;
  }
  return {zmin, zmax};
}
