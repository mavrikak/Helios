/**
 * @file Bsp3D.cpp
 * @brief Implementation of a minimal 3D BSP tree for triangle meshes.
 *
 * Implementation notes:
 * - The splitting plane tolerance is controlled by @c BSPTHRESH.
 * - Memory ownership: triangles referenced by @c bspPolygon::origin are owned
 * elsewhere; this module only allocates/deallocates @c bspPolygon wrappers
 * created during splitting.
 */

#include "Bsp3D.h"
#include <cmath>

/// Distance tolerance for classifying points against a splitter plane.
#define BSPTHRESH 1e-6


// ----------------------------------
// Construction & destruction
// ----------------------------------

Bsp3D::Bsp3D() : bud(false) { elements.clear(); }

Bsp3D::~Bsp3D() { delete bspInside; delete bspOutside; }

// ----------------------------------
// Build operations
// ----------------------------------

// Choose the last pushed element as splitter and remove it from the queue.
bool Bsp3D::SetSplitter() {
  splitter = elements.back();
  elements.pop_back();
  return true;
}

// Wrap a triangle into a working polygon (copy the 3 vertices).
bool Bsp3D::AddElement(Triangle *inElement, int inSide) {
  bspPolygon *poly = new bspPolygon;
  poly->origin = inElement;
  poly->side = inSide;
  for (int i(0); i < 3; i++) poly->nodes.push_back(inElement->Node(i));
  elements.push_back(poly);
  return true;
}

// Add polygon to list
bool Bsp3D::AddElement(bspPolygon *poly) {
  elements.push_back(poly);
  return true;
}

// Set this branch to a bud (leaf): branch with no elements
bool Bsp3D::SetBud() {
  bud = true;
  return true;
}

// Query if this branch is a bud
bool Bsp3D::IsBud() { return bud; }

// ----------------------------------
// Geometry: plane classification
// ----------------------------------

// Test if a point is outside, inside or in plane with this branch's splitter
int Bsp3D::LocTest(rvec test) {
  double dist = LocDist(test);
  if (dist > BSPTHRESH)
    return 1;   // Outside
  else if (dist < -BSPTHRESH)
    return -1;  // Inside
  else
    return 0;   // In plane
}

// Distance to splitter plane
double Bsp3D::LocDist(rvec test) {
  rvec ref = splitter->origin->Center();
  rvec normal = splitter->origin->FBNormalVector(1 - splitter->side);
  rvec diff(test - ref);
  double proj(dot(diff, normal));
  return proj;
}

// ----------------------------------
// Core: split & dispatch
// ----------------------------------

/* Assign poly to proper branch, split if necessary.
 * Returns true if poly has been split into new polygons and can be deleted
 * Returns false if poly remains needed
 */
bool Bsp3D::SplitPoly(bspPolygon *poly) {
  bool inside(false);
  bool outside(false);
  std::vector<int> location;
  for (uint i(0); i < poly->nodes.size(); i++) {
    switch (LocTest(poly->nodes.at(i))) {
      case -1:  // point is inside
        inside = true;
        location.push_back(-1);
        break;
      case 1:   // point is outside
        outside = true;
        location.push_back(1);
        break;
      case 0:   // point is in plane
        location.push_back(0);
    }
  }

  // poly intersects splitter: split into two polygons
  if (inside && outside) {
    // Determine intersection location and add nodes
    if (location.front() * location.back() == -1) {
      double dFirst = LocDist(poly->nodes.front());
      double dLast = LocDist(poly->nodes.back());
      dFirst = dLast / (dFirst - dLast);
      if (dFirst < 0) dFirst = -dFirst;
      dLast = 1. - dFirst;
      poly->nodes.push_back(dFirst * poly->nodes.front() +
                            dLast * poly->nodes.back());
      location.push_back(0);
    }
    for (uint i(0); i < location.size() - 1; i++) {
      if (location.at(i) * location.at(i + 1) == -1) {
        std::vector<rvec>::iterator nodeIter = poly->nodes.begin();
        std::vector<int>::iterator locIter = location.begin();
        for (uint j(0); j < i + 1; j++) {
          nodeIter++;
          locIter++;
        }
        double dFirst = LocDist(poly->nodes.at(i));
        double dLast = LocDist(poly->nodes.at(i + 1));
        dFirst = dLast / (dFirst - dLast);
        if (dFirst < 0) dFirst = -dFirst;
        dLast = 1. - dFirst;
        poly->nodes.insert(nodeIter, 1, dFirst * poly->nodes.at(i) +
                                            dLast * poly->nodes.at(i + 1));
        location.insert(locIter, 1, 0);
      }
    }
    // Create new polygons
    bspPolygon *polyIn = new bspPolygon;
    bspPolygon *polyOut = new bspPolygon;
    polyIn->origin = poly->origin;
    polyIn->side = poly->side;
    polyOut->origin = poly->origin;
    polyOut->side = poly->side;
    for (uint i(0); i < location.size(); i++) {
      switch (location.at(i)) {
        case -1:
          polyIn->nodes.push_back(poly->nodes.at(i));
          break;
        case 1:
          polyOut->nodes.push_back(poly->nodes.at(i));
          break;
        case 0:
          polyIn->nodes.push_back(poly->nodes.at(i));
          polyOut->nodes.push_back(poly->nodes.at(i));
      }
    }
    bspInside->AddElement(polyIn);
    bspOutside->AddElement(polyOut);
    // poly is no longer needed, can be deleted (true)
    return true;
  } else if (inside)
    bspInside->AddElement(poly);
  else if (outside)
    bspOutside->AddElement(poly);
  // poly still needed, do not delete!
  return false;
}

// Perform branching: assign all elements to outside or inside subbranches (or both)
bool Bsp3D::Branch() {
  bspInside = new Bsp3D;
  bspOutside = new Bsp3D;

  for (std::vector<bspPolygon *>::iterator elem = elements.begin();
       elem != elements.end(); elem++) {
    if (SplitPoly(*elem))  // true if poly has been cloned
      delete *elem;
  }

  if ((int)(bspInside->elements.size()) == 0) {
    bspInside->SetBud();
  } else {
    bspInside->SetSplitter();
    bspInside->Branch();
  }
  // ... and the same for the outside branch
  if ((int)(bspOutside->elements.size()) == 0) {
    bspOutside->SetBud();
  } else {
    bspOutside->SetSplitter();
    bspOutside->Branch();
  }
  return true;
}

// ----------------------------------
// Queries
// ----------------------------------

/**
 * @brief Write a node's coordinates to standard output.
 * @param in Input vector containing the node coordinates (x, y, z).
 * @return Always returns false (0).
 * @details Prints the three components of @p in separated by spaces.
 */
bool WriteNode(rvec in) {
  std::cout << in(0) << " " << in(1) << " " << in(2) << " ";
  return 0;
}

// Test if point is inside or outside of branch
int Bsp3D::TestPoint(rvec point) {
  if (LocTest(point) == 1) {  // outside
    if (bspOutside->IsBud()) {
      return 1;
    } else {
      return bspOutside->TestPoint(point);
    }
  } else {  // LocTest(point) <= 0 -> inside (by convention)
    if (bspInside->IsBud()) {
      return -1;
    } else {
      return bspInside->TestPoint(point);
    }
  }
}

