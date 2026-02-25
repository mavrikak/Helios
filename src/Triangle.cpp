/**
 * @file Triangle.cpp
 * @brief Implementation of Triangle geometry caching and utilities.
 */

#include "Triangle.h"
#include <iostream>

// ---------------- Constructors ----------------
// Construct with nodes and domain indices; compute geometry cache.
Triangle::Triangle(rvec *n1, rvec *n2, rvec *n3, int dIndex1, int dIndex2) {
  nodes.push_back(n1);
  nodes.push_back(n2);
  nodes.push_back(n3);

  domainIndeces.push_back(dIndex1);
  domainIndeces.push_back(dIndex2);

  CalculateGeometryData();
}

// Construct with nodes; domains can be set later.
Triangle::Triangle(rvec *n1, rvec *n2, rvec *n3) {
  nodes.push_back(n1);
  nodes.push_back(n2);
  nodes.push_back(n3);

  CalculateGeometryData();
}

// Set front/back domain indices; recompute geometry cache.
int Triangle::SetDomainIndeces(int dIndex1, int dIndex2) {
  domainIndeces.push_back(dIndex1);
  domainIndeces.push_back(dIndex2);
  CalculateGeometryData();
  return 0;
}

// Swap first two nodes (invert orientation).
int Triangle::Invert() {
  rvec *nTmp;
  nTmp = nodes.at(0);
  nodes.at(0) = nodes.at(1);
  nodes.at(1) = nTmp;
  return 0;
}

// ---------------- Geometry cache ----------------
/**
 * @details Recompute normal (FRONT), area, perimeter, centroid, and I1/I2s/I2d.
 *
 * Normal and area:
 *   @f$ \mathbf{n} = (\mathbf{n_2}-\mathbf{n_1}) \times (\mathbf{n_3}-\mathbf{n_1}), 
 *   A = |\mathbf{n}|/2, \hat{\mathbf{n}} = \mathbf{n}/|\mathbf{n}| @f$.
 *
 * Perimeter and centroid:
 *   @f$ \mathbf{r_{cnt}} = (\mathbf{n_1}+\mathbf{n_2}+\mathbf{n_3})/3, 
 *       \Pi = \sum_{i=1}^3 |\mathbf{e_i}|@f$ (sum of edge lengths).
 *
 * Analytic integrals (I1,I2s,I2d) follow the project's closed forms.
 */
int Triangle::CalculateGeometryData() {
  normalVector = cross((rvec)(*(nodes[1]) - *(nodes[0])),
                       (rvec)(*(nodes[2]) - *(nodes[0])));
  surfaceArea = sqrt(dot(normalVector, normalVector));
  normalVector /= surfaceArea;
  surfaceArea /= 2.;
  center = 0;
  perimeter = 0;
  double l[3];
  for (int i(0); i <= 2; ++i) {
    center += Node(i) / 3.;
    l[i] = sqrt(dot(Node((i + 1) % 3) - Node((i + 2) % 3),
                    Node((i + 1) % 3) - Node((i + 2) % 3)));
    perimeter += l[i];
  }
  I1 = 0;
  double &A = surfaceArea, p = perimeter / 2;
  for (int i(0); i <= 2; ++i) {
    double a = l[i], b = l[(i + 1) % 3], c = l[(i + 2) % 3];
    I1 += log(1 - a / p) / a;
    I2s[i] =
        A * A / 30 *
        ((10 + 3 * (c * c - a * a) / (b * b) + 3 * (b * b - a * a) / (c * c)) *
             a -
         (5 - 3 * (a * a - b * b) / (c * c) - 2 * (b * b - c * c) / (a * a)) *
             b -
         (5 - 3 * (a * a - c * c) / (b * b) - 2 * (c * c - b * b) / (a * a)) *
             c +
         (a * a - 3 * b * b - 3 * c * c - 8 * A * A / (a * a)) * 2 *
             log(1 - a / p) / a +
         (a * a - 2 * b * b - 4 * c * c + 6 * A * A / (b * b)) * 4 *
             log(1 - b / p) / b +
         (a * a - 2 * c * c - 4 * b * b + 6 * A * A / (c * c)) * 4 *
             log(1 - c / p) / c);
    I2d[i] =
        A * A / 60 *
        ((-10 + (c * c - a * a) / (b * b) + (b * b - a * a) / (c * c)) * a +
         (5 + (a * a - b * b) / (c * c) - 6 * (b * b - c * c) / (a * a)) * b +
         (5 + (a * a - c * c) / (b * b) - 6 * (c * c - b * b) / (a * a)) * c +
         (2 * a * a - b * b - c * c + 4 * A * A / (a * a)) * 12 *
             log(1 - a / p) / a +
         (9 * a * a - 3 * b * b - c * c + 4 * A * A / (b * b)) * 2 *
             log(1 - b / p) / b +
         (9 * a * a - 3 * c * c - b * b + 4 * A * A / (c * c)) * 2 *
             log(1 - c / p) / c);
  }
  I1 *= -4. / 3 * A * A;
  return 0;
}

// ---------------- Cached getters & simple ops ----------------
double Triangle::getI1() { return I1; }

double Triangle::getI2(rvec *r) {
  for (int i = 0; i < 3; i++) {
    if (nodes[i] == r) return I2s[i];
  }
  exit(1);
  return 0;
}

double Triangle::getI2(rvec *r1, rvec *r2) {
  for (int i = 0; i < 3; i++) {
    if (nodes[i] != r1 and nodes[i] != r2) return I2d[i];
  }
  exit(1);
  return 0;
}

/**
 * @details I2 for a translated configuration (periodic images).
 *
 * If r1 translated by t coincides with r2, returns the single-vertex value I2s;
 * otherwise returns the double-vertex value I2d associated with the remaining node.
 */
double Triangle::getI2Translated(rvec r1, rvec r2, rvec t) {
  rvec r1t(r1 - t);
  const double eps(1e-20);
  if (dot(r1t - r2, r1t - r2) < eps) {
    for (int i = 0; i < 3; i++) {
      if (dot(*nodes[i] - r2, *nodes[i] - r2) < eps) return I2s[i];
    }
  } else {
    for (int i = 0; i < 3; i++) {
      if (dot(*nodes[i] - r1t, *nodes[i] - r1t) > eps and
          dot(*nodes[i] - r2, *nodes[i] - r2) > eps)
        return I2d[i];
    }
  }
#ifdef DEBUG
  std::cout << r1 << r2 << t << "\n";
  for (int i = 0; i < 3; i++) {
    std::cout << *nodes[i] << "\n";
  }
#endif
  exit(1);
  return 0;
}

// ---------------- Basic accessors ----------------
// Node value by index.
rvec Triangle::Node(int index) { return *(nodes[index]); }

// Node pointer by index.
rvec *Triangle::NodePtr(int index) { return nodes[index]; }

// Return FRONT/BACK if dIndex matches, otherwise −1.
int Triangle::DomainPosition(int dIndex) {
  if (domainIndeces[FRONT] == dIndex)
    return FRONT;
  else if (domainIndeces[BACK] == dIndex)
    return BACK;
  else
    return -1;
}

// Reference unit normal (points to FRONT).
rvec Triangle::RefNormalVector() { return normalVector; }

// FRONT/BACK normal.
rvec Triangle::FBNormalVector(int side) {
  if (side == FRONT)
    return normalVector;
  else
    return -normalVector;
}

// Area (cached).
double Triangle::Area() { return surfaceArea; }

// Perimeter (cached).
double Triangle::Perimeter() { return perimeter; }

// Centroid (cached).
rvec Triangle::Center() { return center; }

// Pair {front, back} domain indices.
std::pair<int, int> Triangle::BorderingIndeces() {
  return std::pair<int, int>(domainIndeces.at(0), domainIndeces.at(1));
}

// Replace node pointers and recompute geometry cache.
int Triangle::AssignNodes(rvec *n1, rvec *n2, rvec *n3) {
  nodes.clear();
  nodes.push_back(n1);
  nodes.push_back(n2);
  nodes.push_back(n3);
  CalculateGeometryData();
  return 0;
}

// Analytic R^2 integral over triangle (self term), used in some kernels.
double Triangle::IntRR() {
  rvec xp, yp, zp, rc;
  double jacobian;
  double mA[3][3], mB[3][3];
  double result;

  xp = *(nodes.at(1)) - *(nodes.at(0));
  yp = *(nodes.at(2)) - *(nodes.at(0));
  zp = cross(xp, yp);
  jacobian = sqrt(dot(zp, zp));
  for (int i(0); i < 3; i++) {
    mA[i][0] = xp(i);
    mA[i][1] = yp(i);
    mA[i][2] = zp(i);
  }
  for (int i(0); i < 3; i++) {
    for (int j(0); j < 3; j++) {
      mB[i][j] = 0.;
      {
        for (int k(0); k < 3; k++) mB[i][j] += mA[k][i] * mA[k][j];
      }
    }
  }
  rc = rvec((mA[0][0] + mA[0][1]) / 3., (mA[1][0] + mA[1][1]) / 3.,
            (mA[2][0] + mA[2][1]) / 3.);
  result = jacobian *
           ((mB[0][0] + mB[1][1]) / 12. + (mB[0][1] + mB[1][0]) / 24. +
            dot(*(nodes.at(0)), rc) + dot(*(nodes.at(0)), *(nodes.at(0))) / 2.);
  return result;
}

// Replace domain index n1 with n2 if present.
int Triangle::SwitchIndex(int n1, int n2) {
  if (domainIndeces[0] == n1) domainIndeces[0] = n2;
  if (domainIndeces[1] == n1) domainIndeces[1] = n2;
  return 0;
}

// ---------------- Adjacency / common entities ----------------
// Count shared nodes with triangle T (0..3).
int Triangle::FindAdjacency(Triangle *T) {
  int ret = 0;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (NodePtr(i) == T->NodePtr(j)) {
        ret++;
        break;
      }
    }
  }
  return ret;
}

// Find common edge (local indices) with T; return {i+1, -(j+1)} or {0,0}.
std::pair<int,int> Triangle::FindCommonEdgeIndices(Triangle* T) {
  for (int i = 0; i < 3; ++i) {
    int iEdgeNode1 = i, iEdgeNode2 = (i + 1) % 3;
    for (int j = 0; j < 3; ++j) {
      int jEdgeNode1 = j, jEdgeNode2 = (j + 1) % 3;
      if ( ( this->NodePtr(iEdgeNode1) == T->NodePtr(jEdgeNode1) && 
             this->NodePtr(iEdgeNode2) == T->NodePtr(jEdgeNode2) ) ||
           ( this->NodePtr(iEdgeNode1) == T->NodePtr(jEdgeNode2) && 
             this->NodePtr(iEdgeNode2) == T->NodePtr(jEdgeNode1) ) ) {
        return { i + 1, -(j + 1) };
      }
    }
  }
  return { 0, 0 }; // no common edge
}

// Find common vertex (local indices) with T; return {i+1, j+1} or {0,0}.
std::pair<int,int> Triangle::FindCommonVertexIndices(Triangle* T) {
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      if (this->NodePtr(i) == T->NodePtr(j)) {
        return { i + 1, j + 1 };
      }
    }
  }
  return { 0, 0 }; // no common edge
}

// Count shared nodes with T after translating T by vector t.
int Triangle::FindAdjacency(Triangle *T, rvec t) {
  int ret = 0;
  double eps = 1e-4;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      rvec rij(Node(i) - t - T->Node(j));
      if (dot(rij, rij) < eps) {
        ret++;
        break;
      }
    }
  }
  return ret;
}

// ---------------- Point-in-triangle ----------------
/**
 * @details Check if a point lies on/near the triangle (plane-project + barycentric).
 *
 * Steps:
 *  1) Reject if @f$ |(\mathbf{p} − \mathbf{n_0}) \cdot \mathbf{n}| > tol @f$ (off the plane).
 *  2) Build an orthonormal basis {unitX, unitY} in the triangle plane.
 *  3) Compute simplex (barycentric) coordinates (zeta) in 2D.
 *  4) Accept if zeta_i >= −tol and sum(zeta) <= 1+tol.
 */
bool Triangle::containsPoint(const rvec &p, double tol)
{
  double A = Area();

  // Degenerate triangle
  if (std::fabs(A) < tol) return false;
  
  // Coplanarity: distance from p to the triangle's plane
  double dist = dot(p - Node(0), normalVector);
  if (std::fabs(dist) > tol) return false;

  // unit vectors
  rvec unitX = Node(1) - Node(0);  unitX /= sqrt( dot(unitX, unitX) );
  rvec temp  = Node(2) - Node(0);
  rvec unitZ = cross(unitX, temp); unitZ /= sqrt( dot(unitZ, unitZ) );
  rvec unitY = cross(unitZ, unitX);

  // vertices in triangle plane
  std::vector<rvec> nodesInPlane;
  nodesInPlane.push_back( rvec( 0.0, 0.0, 0.0 ) );                  // Node(0) in plane
  temp = Node(1) - Node(0);
  nodesInPlane.push_back( rvec( clean( dot(temp, unitX) ), 
                                clean( dot(temp, unitY) ), 0.0 ) ); // Node(1) in plane
  temp = Node(2) - Node(0);
  nodesInPlane.push_back( rvec( clean( dot(temp, unitX) ), 
                                clean( dot(temp, unitY) ), 0.0 ) ); // Node(2) in plane

  // Observation point in triangle plane
  rvec pos( clean( dot( unitX, p - Node(0) ) ), 
            clean( dot( unitY, p - Node(0)) ), 0.0 );

  // Simplex coordinates
  std::vector<double> zeta(3);
  for (int i = 0; i < 3; i++) {
    zeta[i] = 
    ( ( nodesInPlane[(i + 1) % 3][0] * nodesInPlane[(i + 2) % 3][1] - 
        nodesInPlane[(i + 2) % 3][0] * nodesInPlane[(i + 1) % 3][1] ) + 
      pos[0] * ( nodesInPlane[(i + 1) % 3][1] - nodesInPlane[(i + 2) % 3][1] ) +
      pos[1] * ( nodesInPlane[(i + 2) % 3][0] - nodesInPlane[(i + 1) % 3][0] ) ) /
    ( 2. * A );
  }

  // Check if point is in triangle (allow a tiny epsilon outside)
  if (zeta[0] < -tol || zeta[1] < -tol ||
      zeta[2] < -tol || (zeta[0] + zeta[1] + zeta[2]) > 1.0 + tol) 
    return false;

  return true;
}