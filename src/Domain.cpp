/**
 * @file Domain.cpp
 * @brief Implementation of the Domain base class.
 */

#include "Domain.h"
#include <iostream>

// -----------------------------
// Construction & lifetime
// -----------------------------
Domain::Domain(SurfaceMesh *mesh, int dIndex, dcmplx inVacuumWavelength)
    : vacuumWavelength(inVacuumWavelength),
      omega(2 * PI / inVacuumWavelength * CVAC),
      index(dIndex) { AddTrianglesFromMesh(mesh, dIndex); }

Domain::~Domain() {}

int Domain::Index() { return index; }

// -----------------------------
// Mesh / basis assembly
// -----------------------------
int Domain::AddTrianglesFromMesh(SurfaceMesh *mesh, int dIndex) {
  std::vector<Triangle *> inTriangles = mesh->BorderingTriangles(dIndex);
  std::vector<Triangle *>::iterator triIter;
  for (triIter = inTriangles.begin(); triIter != inTriangles.end(); triIter++) {
    int side((*triIter)->DomainPosition(dIndex));
    if (side == FRONT || side == BACK) {
      // One RWG per directed edge (three hRWGs per triangle here)
      rwgFuns.push_back(RWGFun(mesh, *triIter, (*triIter)->NodePtr(0), side));
      rwgFuns.push_back(RWGFun(mesh, *triIter, (*triIter)->NodePtr(1), side));
      rwgFuns.push_back(RWGFun(mesh, *triIter, (*triIter)->NodePtr(2), side));

      // Feed triangle to BSP for inside/outside classification
      bspTree.AddElement(*triIter, side);
    }
  }

  // Build the BSP tree (choose last polygon as splitter and branch recursively)
  bspTree.SetSplitter();
  bspTree.Branch();

  std::cout << "Added " << inTriangles.size() << " triangles and "
            << rwgFuns.size() << " RWG functions to domain " << dIndex
            << std::endl;

  return inTriangles.size();
}

// -----------------------------
// Global accuracy flag
// -----------------------------
bool Domain::NeedsAccurate = false;
void Domain::EnableAccurate() { NeedsAccurate = true; }

// -----------------------------
// Formulation assembly hooks
// -----------------------------
int Domain::FillLSEMatrix(blitz::Array<dcmplx, 2> *matrix) {
  return formulation->FillLSEMatrix(matrix);
}

int Domain::FillLSEMatrixAndGradient(blitz::Array<dcmplx, 2> *matrix,
                          blitz::Array<dcmplx, 2> *gradmatrix[3]) {
  return formulation->FillLSEMatrixAndGradient(matrix, gradmatrix);
}

int Domain::FillLSEVector(blitz::Array<dcmplx, 1> *vector) {
  return formulation->FillLSEVector(vector);
}

// -----------------------------
// Geometry queries
// -----------------------------
bool Domain::IsInside(rvec point) {
  if (bspTree.TestPoint(point) == -1) return true;
  return false;
}

// -----------------------------
// Incident field handling
// -----------------------------
bool Domain::hasDipoles() {
  for (std::vector<IncidentField *>::iterator incIter = incidentFields.begin();
        incIter != incidentFields.end(); incIter++) {
    if ((*incIter)->IsDipole()) return true;
  }
  return false;
}

// Add an incident condition
int Domain::AddIncident(IncidentField *inIncident) {
  incidentFields.push_back(inIncident);
  return 0;
}

// -----------------------------
// Field evaluation (post‑processing)
// -----------------------------
// Secondary (scattered) electric field at a point from the solution vector.
cvec Domain::SecFieldE(blitz::Array<dcmplx, 1> *solVector, rvec pos) {
  cvec result(0, 0, 0);
  int offset(solVector->shape()(0) / 2);
  using std::vector;
  using std::pair;
  vector<RWGFun>::iterator rwgM, rwgN;

  // Group hRWGs that belong to the same triangle
  vector<RWGFun *> rwgPerTri;
  for (rwgM = rwgFuns.begin(); rwgM != rwgFuns.end(); rwgM = rwgN) {
    rwgPerTri.clear();
    for (rwgN = rwgM;
         rwgN != rwgFuns.end() and rwgN->TrianglePtr() == rwgM->TrianglePtr();
         rwgN++) {
      rwgPerTri.push_back(&(*rwgN));
    }
    // Compute operators for this group of hRWGs
    vector<pair<cvec, cvec> > DeltaKappa(
        rwgPerTri.size(), make_pair(cvec(0., 0., 0.), cvec(0., 0., 0.)));
    grnFun->SameTriDeltaKappa(pos, rwgPerTri, DeltaKappa);

    // Accumulate field contributions from this group of hRWGs
    for (int i = 0; i < (int)rwgPerTri.size(); ++i) {
      result -= I * omega * MU0 * Mu(pos) * DeltaKappa[i].first *
                    (*solVector)(rwgPerTri[i]->EdgeIndex()) +
                DeltaKappa[i].second *
                    (*solVector)(rwgPerTri[i]->EdgeIndex() + offset);
    }
  }
  return result;
}

// Total incident electric field at a point from all incident sources.
cvec Domain::IncFieldE(rvec pos) {
  cvec incField(0, 0, 0);
  for (std::vector<IncidentField *>::iterator incIter = incidentFields.begin();
       incIter != incidentFields.end(); incIter++)
    incField += (*incIter)->EvaluateE(pos);
  return incField;
}

// Secondary (scattered) magnetic field at a point from the solution vector.
cvec Domain::SecFieldH(blitz::Array<dcmplx, 1> *solVector, rvec pos) {
  cvec result(0, 0, 0);
  int offset(solVector->shape()(0) / 2);
  using std::vector;
  using std::pair;
  vector<RWGFun>::iterator rwgM, rwgN;

  // Group hRWGs that belong to the same triangle
  vector<RWGFun *> rwgPerTri;
  for (rwgM = rwgFuns.begin(); rwgM != rwgFuns.end(); rwgM = rwgN) {
    rwgPerTri.clear();
    for (rwgN = rwgM;
         rwgN != rwgFuns.end() and rwgN->TrianglePtr() == rwgM->TrianglePtr();
         rwgN++) {
      rwgPerTri.push_back(&(*rwgN));
    }
    // Compute operators for this group of hRWGs
    vector<pair<cvec, cvec> > DeltaKappa(
        rwgPerTri.size(), make_pair(cvec(0., 0., 0.), cvec(0., 0., 0.)));
    grnFun->SameTriDeltaKappa(pos, rwgPerTri, DeltaKappa);

    // Accumulate field contributions from this group of hRWGs
    for (int i = 0; i < (int)rwgPerTri.size(); ++i) {
      result -= I * omega * EPS0 * Epsilon(pos) * DeltaKappa[i].first *
                    (*solVector)(rwgPerTri[i]->EdgeIndex() + offset) -
                DeltaKappa[i].second * (*solVector)(rwgPerTri[i]->EdgeIndex());
    }
  }
  return result;
}

// Secondary E-/H-fields at a point.
std::pair<cvec, cvec> Domain::SecFieldEandH(blitz::Array<dcmplx, 1> *solVector,
                                            rvec pos) {
  std::pair<cvec, cvec> result(make_pair(cvec(0, 0, 0), cvec(0, 0, 0)));
  int offset(solVector->shape()(0) / 2);
  using std::vector;
  using std::pair;
  vector<RWGFun>::iterator rwgM, rwgN;

  // Group hRWGs that belong to the same triangle
  vector<RWGFun *> rwgPerTri;
  for (rwgM = rwgFuns.begin(); rwgM != rwgFuns.end(); rwgM = rwgN) {
    rwgPerTri.clear();
    for (rwgN = rwgM;
         rwgN != rwgFuns.end() and rwgN->TrianglePtr() == rwgM->TrianglePtr();
         rwgN++) {
      rwgPerTri.push_back(&(*rwgN));
    }
    // Compute operators for this group of hRWGs
    vector<pair<cvec, cvec> > DeltaKappa(
        rwgPerTri.size(), make_pair(cvec(0., 0., 0.), cvec(0., 0., 0.)));
    grnFun->SameTriDeltaKappa(pos, rwgPerTri, DeltaKappa);

    // Accumulate field contributions from this group of hRWGs
    for (int i = 0; i < (int)rwgPerTri.size(); ++i) {
      result.first -= I * omega * MU0 * Mu(pos) * DeltaKappa[i].first *
                          (*solVector)(rwgPerTri[i]->EdgeIndex()) +
                      DeltaKappa[i].second *
                          (*solVector)(rwgPerTri[i]->EdgeIndex() + offset);
      result.second -=
          I * omega * EPS0 * Epsilon(pos) * DeltaKappa[i].first *
              (*solVector)(rwgPerTri[i]->EdgeIndex() + offset) -
          DeltaKappa[i].second * (*solVector)(rwgPerTri[i]->EdgeIndex());
    }
  }
  return result;
}

// Total incident magnetic field at a point from all incident sources.
cvec Domain::IncFieldH(rvec pos) {
  cvec incField(0, 0, 0);
  for (std::vector<IncidentField *>::iterator incIter = incidentFields.begin();
       incIter != incidentFields.end(); incIter++)
    incField += (*incIter)->EvaluateH(pos);
  return incField;
}

// -----------------------------
// Power computation
// -----------------------------
/**
 * @brief Compute the element-wise complex conjugate of a 3-vector.
 * @param in Input complex vector.
 * @return Complex vector containing conjugated components of @p in.
 */
cvec conj(cvec in) { return cvec(conj(in(0)), conj(in(1)), conj(in(2))); }

/**
 * @brief Extract the real part of each component of a complex 3-vector.
 * @param in Input complex vector.
 * @return Real vector with the real parts of @p in.
 */
rvec real(cvec in) { return rvec(real(in(0)), real(in(1)), real(in(2))); }

/* Calculation details based on Section 4 of 
 * A. Kern, "Realistic modeling of 3D plasmonic
 * systems: A surface integral equation approach,"
 * Ph.D. thesis, Nanophotonics and Metrology Lab.,
 * EPFL, Lausanne, Switzerland, 2011.
 */
double Domain::Power(blitz::Array<dcmplx, 1> *solVector) {
  double result(0.);
  std::vector<Triangle *> vTri;
  std::vector<dcmplx> vA, vB;
  std::vector<cvec> vAP, vBP;
  std::vector<rvec> vNormals;

  for (std::vector<RWGFun>::iterator rwgIter = rwgFuns.begin();
       rwgIter != rwgFuns.end(); rwgIter++) {
    int triIndex(-1);
    for (uint i(0); i < vTri.size(); i++)
      if (vTri.at(i) == rwgIter->TrianglePtr()) triIndex = i;
    if (triIndex == -1) {
      triIndex = vTri.size();
      vTri.push_back(rwgIter->TrianglePtr());
      vA.push_back(0.);
      vB.push_back(0.);
      vAP.push_back(cvec(0, 0, 0));
      vBP.push_back(cvec(0, 0, 0));
      vNormals.push_back(
          rwgIter->TrianglePtr()->FBNormalVector(rwgIter->Side()));
    }
    int offset(solVector->shape()(0) / 2);
    dcmplx addA(conj((*solVector)(rwgIter->EdgeIndex()) *
                     rwgIter->RWGPreFactor() / 2.));
    dcmplx addB((*solVector)(rwgIter->EdgeIndex() + offset) *
                rwgIter->RWGPreFactor() / 2.);
    vA.at(triIndex) = vA.at(triIndex) + addA;
    vB.at(triIndex) = vB.at(triIndex) + addB;
    vAP.at(triIndex) = vAP.at(triIndex) + addA * rwgIter->FreeVertex();
    vBP.at(triIndex) = vBP.at(triIndex) + addB * rwgIter->FreeVertex();
  }
  for (uint i(0); i < vTri.size(); i++) {
    rvec normal = vNormals.at(i);
    rvec centroid = vTri.at(i)->Center();
    double area = vTri.at(i)->Area();
    result +=
        real(area *
             dot(normal, cross((cvec)centroid, (cvec)(vB.at(i) * vAP.at(i) -
                                                      vA.at(i) * vBP.at(i))) +
                             cross((cvec)vAP.at(i), (cvec)vBP.at(i)))) /
        2.;
  }
  return result;
}

// -----------------------------
// Incident cleanup
// -----------------------------
int Domain::ClearIncidentFields() {
  while (incidentFields.size() > 0) {
    delete incidentFields.back();
    incidentFields.pop_back();
  }
  return 0;
}

// -----------------------------
// Debug utilities
// -----------------------------
int Domain::PrintAllRWGCentroids() {
  std::vector<RWGFun>::iterator rwgIter;
  for (rwgIter = rwgFuns.begin(); rwgIter != rwgFuns.end(); rwgIter++) {
    Triangle *ptr = rwgIter->TrianglePtr();
    rvec n1 = ptr->Node(0);
    rvec n2 = ptr->Node(1);
    rvec n3 = ptr->Node(2);
    rvec ctr = (n1 + n2 + n3) / 3.;
    std::cout << "Triangle: " << n1 << std::endl;
    std::cout << "          " << n2 << std::endl;
    std::cout << "          " << n3 << std::endl;
    std::cout << "     Fvx: " << rwgIter->FreeVertex() << std::endl;
    std::cout << "     RWG: " << rwgIter->Evaluate(ctr) << std::endl;
  }
  return 0;
}

int Domain::PrintAllEdges(SurfaceMesh *mesh) {
  for (int i(0); i < mesh->EdgeCount(); ++i) {
    std::cout << "EDGE " << i << " :" << std::endl;
    for (unsigned int j(0); j < rwgFuns.size(); ++j) {
      if (rwgFuns[j].EdgeIndex() == i) {
        Triangle *ptr = rwgFuns[j].TrianglePtr();
        rvec n1 = ptr->Node(0);
        rvec n2 = ptr->Node(1);
        rvec n3 = ptr->Node(2);
        rvec ctr = (n1 + n2 + n3) / 3.;

        std::cout << "  Center: " << ctr << std::endl;
        std::cout << "     RWG: " << rwgFuns[j].Evaluate(ctr) << std::endl;
      }
    }
  }
  return 0;
}

RWGFun *Domain::RGWFunPtr(int index) { return &(rwgFuns[index]); }

int Domain::InitLookupTable(LookupTableBin *tin) {
  grnFun->InitLookupTable(tin);
  return 0;
}
