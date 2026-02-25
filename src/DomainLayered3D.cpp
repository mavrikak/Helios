/**
 * @file DomainLayered3D.cpp
 * @brief Implementation of the layered 3D domain.
 */

#include "DomainLayered3D.h"
#include <iostream>
#include "Dipole.h"
#include "IncidentField.h"
#include "PlaneWave.h"
#include "SIEFormPMCHW.h"
#include "SommerfeldIntegrator.h"
#include "iofunctions.h"
#include "Gaussian.h"
#include <limits>



// -----------------------------------------------------------------------------
// Helper: read XYZ positions from a text file (one triplet per line)
// -----------------------------------------------------------------------------
/**
 * @brief Reads a set of 3D positions from a file.
 *
 * This function reads a file containing 3D coordinates (x, y, z) line by line
 * and returns a vector of rvec objects representing these positions.
 *
 * @param filename The name of the file containing position data.
 * @return A vector of rvec objects containing the parsed positions.
 * @throws std::runtime_error if the file cannot be opened.
 */
std::vector<rvec> readPositionsFromFile(const std::string& filename) {
  std::vector<rvec> positions;
  std::ifstream file(filename);

  if (!file.is_open()) {
      throw std::runtime_error("Could not open file: " + filename);
  }

  std::string line;
  while (std::getline(file, line)) {
      std::istringstream iss(line);
      double x, y, z;
      if (iss >> x >> y >> z) {
          positions.emplace_back(x, y, z);
      }
  }

  file.close();
  return positions;
}

// -----------------------------------------------------------------------------
// Construction & lifetime
// -----------------------------------------------------------------------------
DomainLayered3D::DomainLayered3D(SurfaceMesh *inMesh, int inDIndex, 
                                 const std::vector<dcmplx>& inEpsilon, 
                                 const std::vector<dcmplx>& inMu, 
                                 const std::vector<double>& inZVals,
                                 dcmplx inVacuumWavelength)
    : Domain(inMesh, inDIndex, inVacuumWavelength),
      epsilon(inEpsilon),
      mu(inMu),
      zValsInterfaces(inZVals),
      grnFunLayered(nullptr),
      grnFunLayeredFields(nullptr) {
  // Derive layer thicknesses from interface positions
  thickness.resize(zValsInterfaces.size() - 1);
  for (size_t i = 1; i < zValsInterfaces.size(); ++i) {
      thickness[i - 1] = zValsInterfaces[i] - zValsInterfaces[i - 1];
  }

  // Collect triangle sample points for grid sizing (Dunavant rule per triangle)
  for (size_t i = 0; i < rwgFuns.size(); i++) {
    rvec r[3];
    for (int j = 0; j < 3; j++) { r[j] = rwgFuns[i].TrianglePtr()->Node(j); }
    int Ndunavant; double(*xdunavant)[3];
    Ndunavant = NeedsAccurate ? dunavant::N17 : dunavant::N1;
    xdunavant = NeedsAccurate ? dunavant::x17 : dunavant::x1;
    for (int k = 0; k < Ndunavant; ++k) {
      triPoints.push_back(r[0] * xdunavant[k][0] + r[1] * xdunavant[k][1] +
                          r[2] * xdunavant[k][2]);
    }
    // Add the triangle vertices to triPoints
    //triPoints.push_back(r[0]);
    //triPoints.push_back(r[1]);
    //triPoints.push_back(r[2]);
  }

  // Slice & grid for tabulation
  auto pointGroups = slice(triPoints, &triPoints);  
  tabulationGrids = computeGrid(pointGroups);

  // Decide if line‑integral terms are needed (interfaces cutting the mesh extent)
  enableLineInt = false;
  std::pair<double, double> zMinMax = inMesh->zMinMax();
  for (size_t i = 0; i < zValsInterfaces.size(); ++i) {
    if (zMinMax.first < zValsInterfaces[i] - 1e-10 && 
        zMinMax.second > zValsInterfaces[i] + 1e-10) {
      enableLineInt = true;
      break;
    }
  }

  // Calculate midLayerIndex based on mesh bounds
  double zMid = (zMinMax.first + zMinMax.second) * 0.5;
  // Use logic similar to basic GetLayerIndex but without interface ambiguity for mid-point
  midLayerIndex = -1;
  if (zMid < zValsInterfaces.front()) midLayerIndex = 0;
  else if (zMid > zValsInterfaces.back()) midLayerIndex = zValsInterfaces.size();
  else {
    for(size_t i = 0; i < zValsInterfaces.size() - 1; ++i) {
       if(zMid >= zValsInterfaces[i] && zMid < zValsInterfaces[i + 1]) {
           midLayerIndex = i + 1; break;
       }
    }
  }
  if (midLayerIndex == -1) midLayerIndex = 0; // Fallback  

  // Build Green's functions and helpers
  InitGrnFun();
  k_L = grnFunLayered->GetLayersWavenumbers();
  layeredUtils = new LayeredMediaUtils(k_L, epsilon, mu, zValsInterfaces, thickness, midLayerIndex);
}

DomainLayered3D::DomainLayered3D(SurfaceMesh *inMesh, int inDIndex, 
                                 const std::vector<dcmplx>& inEpsilon, 
                                 const std::vector<dcmplx>& inMu, 
                                 const std::vector<double>& inZVals, 
                                 dcmplx inVacuumWavelength, bool parent)
    : Domain(inMesh, inDIndex, inVacuumWavelength),
      epsilon(inEpsilon),
      mu(inMu),
      zValsInterfaces(inZVals),
      grnFunLayered(nullptr),
      grnFunLayeredFields(nullptr) {
  if (!parent) InitGrnFun();
}

/**
 * @note #grnFun is an array allocated with new[] in InitGrnFun(), so delete[] is used here.
 */
DomainLayered3D::~DomainLayered3D() {
  delete[] grnFun;            // per‑layer homogeneous GFs
  delete grnFunLayered;       // layered GF
  delete grnFunLayeredFields; // layered GF for fields
  delete layeredUtils;        // layered media utilities
}


// -----------------------------------------------------------------------------
// Green's functions and formulation
// -----------------------------------------------------------------------------
int DomainLayered3D::InitGrnFun() {
  // One homogeneous GF per layer (intra‑layer pieces)
  grnFun = new GreenFHom3D[zValsInterfaces.size() + 1];
  for (size_t i = 0; i < zValsInterfaces.size() + 1; ++i) {
      grnFun[i] = GreenFHom3D(vacuumWavelength, epsilon[i], mu[i]);
  }

  // Layered GF (Sommerfeld integral with tabulation)
  grnFunLayered = new GreenFLayered3D(vacuumWavelength, epsilon, mu,
                                      zValsInterfaces, thickness, 
                                      tabulationGrids, enableLineInt,
                                      midLayerIndex);
  return 0;
}

void DomainLayered3D::newGrnFunLayered(std::vector<rvec> &posvec) {
  // Keep only observation points inside this domain
  std::vector<rvec> obsPoints;
  for (size_t i = 0; i < posvec.size(); i++) {
    if ( this->IsInside(posvec[i]) ) obsPoints.push_back(posvec[i]);
  }

  // Build grids tailored to these points and create a fresh layered GF for fields
  auto pointGroups = slice(obsPoints, &triPoints);
  tabulationGridsFields = computeGrid(pointGroups);
  grnFunLayeredFields = 
      new GreenFLayered3D(vacuumWavelength, epsilon, mu,
                          zValsInterfaces, thickness, 
                          tabulationGridsFields, false, midLayerIndex);

  // If there are dipoles, wire each to its own evaluation GF with a grid that
  // includes the dipole location (to avoid interpolation artifacts near source)
  for (std::vector<IncidentField *>::iterator incIter = incidentFields.begin();
       incIter != incidentFields.end(); incIter++) {
    if ( (*incIter)->IsDipole() ) {
      Dipole* dipPtr = dynamic_cast<Dipole*>(*incIter);
      if (dipPtr) {
        std::vector<rvec> dipPoint = { dipPtr->Location() };
        auto fieldPointGroups = slice(obsPoints, &dipPoint);
        std::vector<Grid> tabGridFields = computeGrid(fieldPointGroups);
        GreenF* fGrnFields = 
            new GreenFLayered3D(vacuumWavelength, epsilon, mu,
                                zValsInterfaces, thickness, tabGridFields,
                                false, midLayerIndex);
        dipPtr->setGrnFunLayered(fGrnFields, tabGridFields);
      }
      else {
        std::cerr << "Error: Dynamic cast to Dipole failed." << std::endl;
      }
    }
  }
}

int DomainLayered3D::InitFormulation(dcmplx inVacWavelength) {
  formulation = new SIEFormPMCHW(&rwgFuns, grnFun, &incidentFields, inVacWavelength,
                                 grnFunLayered, epsilon, mu, zValsInterfaces, 
                                 midLayerIndex);
  return 0;
}

// -----------------------------------------------------------------------------
// Layer indexing
// -----------------------------------------------------------------------------
int DomainLayered3D::GetLayerIndex(rvec pos) const {
  double z = pos[2];  // Extract the z-component of the position
  int index = -1;
  const double eps = 1e-6; // Small tolerance for floating-point comparisons

  // Check if point is on an interface
  for (size_t i = 0; i < zValsInterfaces.size(); ++i) {
    if (std::abs(z - zValsInterfaces[i]) <= eps) {
      // Interface i separates layer i and layer i+1
      int l_lower = i;
      int l_upper = i + 1;
      index = std::abs(l_lower - midLayerIndex) < std::abs(l_upper - midLayerIndex) 
            ? l_lower : l_upper; break;
    }
  }

  for (size_t i = 0; i < zValsInterfaces.size() - 1; ++i) {
    if (z > zValsInterfaces[i] + eps && z < zValsInterfaces[i + 1] - eps) {
      index = i + 1; break;
    }
  }

  // If z is beyond the last interface, return the index of the last layer
  if (z > zValsInterfaces.back() + eps) {
    index = zValsInterfaces.size();
  }
  
  // If z is below the first interface, return the index of the first layer
  if (z < zValsInterfaces.front() - eps) {
    index = 0;
  }

  return index;
}

// -----------------------------------------------------------------------------
// Layer‑optics helpers (delegations to LayeredMediaUtils)
// -----------------------------------------------------------------------------
std::vector<dcmplx> DomainLayered3D::FresnelCoeff(int layerIndexI, int layerIndexJ, 
                                                  dcmplx krho, const std::string& waves) const {
  return this->layeredUtils->FresnelCoeff(layerIndexI, layerIndexJ, krho, waves);
}

std::vector<std::vector<dcmplx>> DomainLayered3D::InterfaceTransferMatrix(
    int layerIndexI, int layerIndexJ, dcmplx krho, const std::string& waves) const {
  return this->layeredUtils->InterfaceTransferMatrix(layerIndexI, layerIndexJ, krho, waves);
}

/* External LAPACK functions for solving linear systems. These functions 
 * are used for LU factorization (`zgetrf_`) and solving linear 
 * equations using the LU decomposition (`zgetrs_`).
 */
extern "C" {
  /**
  * @brief Computes the LU factorization of a general M-by-N matrix.
  *
  * @param dim1 Number of rows of the matrix (for `zgetrf_`).
  * @param dim2 Number of columns of the matrix (for `zgetrf_`).
  * @param a Pointer to the matrix to be factored or solved.
  * @param lda Leading dimension of the matrix.
  * @param ipiv Pointer to the pivot indices.
  * @param info Output status (0 if successful).
  */
  void zgetrf_(int *dim1, int *dim2, dcmplx *a, int *lda, int *ipiv, int *info);

  /**
  * @brief Solves a system of linear equations using the LU factorization.
  *
  * @param TRANS Specifies the form of the system of equations ('N' for no transpose, etc.).
  * @param N Order of the coefficient matrix.
  * @param NRHS Number of right-hand sides.
  * @param A Pointer to the LU-factored coefficient matrix.
  * @param LDA Leading dimension of A.
  * @param IPIV Pointer to the pivot indices from `zgetrf_`.
  * @param B Pointer to the right-hand side vector(s).
  * @param LDB Leading dimension of B.
  * @param INFO Output status (0 if successful).
  */
  void zgetrs_(char *TRANS, int *N, int *NRHS, dcmplx *A, int *LDA, int *IPIV,
              dcmplx *B, int *LDB, int *INFO);
}

blitz::Array<dcmplx, 2> DomainLayered3D::Secondary(dcmplx krho, int observationLayer, 
                                                   int sourceLayer, const std::string& waves) {
  return this->layeredUtils->Secondary(krho, observationLayer, sourceLayer, waves);
}

std::vector<dcmplx> DomainLayered3D::Fields(dcmplx krho, std::vector<double> observationPointsZ, 
                                            double sourcePointZ, const std::string& waves) {
  return this->layeredUtils->Fields(krho, observationPointsZ, sourcePointZ, waves);
}

int DomainLayered3D::writeFields(dcmplx krho, double start, double end, int numPoints, 
                                 double sourcePointZ, const std::string& waves) {
  std::vector<double> zPoints(numPoints);
  double step = (end - start) / (numPoints - 1);
  for (int i = 0; i < numPoints; ++i) {
      zPoints[i] = start + i * step;
  }

  // Compute the fields for all points
  std::vector<dcmplx> fields = Fields(krho, zPoints, sourcePointZ, waves);

  // Construct file name with waves information
  std::string fileName = "FieldsOutput_" + waves + ".txt";

  // Open a file to write the results
  std::ofstream outFile(fileName);
  if (!outFile) {
      std::cerr << "Error: Could not open file for writing." << std::endl;
      return 1;
  }

  // Write header for MATLAB readability
  outFile << "% z  RealPart  ImagPart\n";

  // Write the Fields results to the file
  for (int i = 0; i < zPoints.size(); ++i) {
      outFile << zPoints[i] << " " << std::real(fields[i]) << " " << std::imag(fields[i]) << "\n";
  }
  outFile.close();
  std::cout << "Fields data written to " << fileName << std::endl;

  return 0;
}

// -----------------------------------------------------------------------------
// Material queries
// -----------------------------------------------------------------------------
dcmplx DomainLayered3D::Epsilon(rvec pos) { return epsilon[GetLayerIndex(pos)]; }

dcmplx DomainLayered3D::Mu(rvec pos) { return mu[GetLayerIndex(pos)]; }

// -----------------------------------------------------------------------------
// Incident fields (layered variants)
// -----------------------------------------------------------------------------
int DomainLayered3D::AddPlaneWave(dcmplx /*wavelength*/, rvec propagationDirection,
                                  cvec polarization) {
  IncidentField *ptr = new PlaneWave(k_L, epsilon, mu, zValsInterfaces, thickness,
                                     propagationDirection, polarization);
  incidentFields.push_back(ptr);
  return 0;
}

int DomainLayered3D::AddDipole(rvec location, cvec polarization) {
  std::vector<rvec> locPoints = { location };
  auto tempPointGroups = slice(triPoints, &locPoints);
  std::vector<Grid> tempTabulationGrids = computeGrid(tempPointGroups);
  GreenF* tempGrnFunHom = 
        new GreenFHom3D(vacuumWavelength, Epsilon(location), Mu(location));
  GreenF* tempGrnFunLayered = 
        new GreenFLayered3D(vacuumWavelength, epsilon, mu,
                            zValsInterfaces, thickness, tempTabulationGrids,
                            false, midLayerIndex);
  LayeredMediaUtils* tempLayeredUtils = 
        new LayeredMediaUtils(k_L, epsilon, mu, zValsInterfaces, 
                              thickness, midLayerIndex);
  IncidentField *ptr =
        new Dipole(location, polarization, vacuumWavelength, tempGrnFunHom, 
                   Mu(location), tempTabulationGrids, tempGrnFunLayered,
                   tempLayeredUtils);
  incidentFields.push_back(ptr);
  return 0;
}

// Not implemented for layered media yet
int DomainLayered3D::AddGaussian(double wavelength, double waist, rvec focal_point,
                                 rvec propagationDirection, cvec polarization) {
  IncidentField *ptr =
      new Gaussian(wavelength, waist, focal_point, real(epsilon[0]), // TEMPORARY CHANGE
                   propagationDirection, polarization);
  incidentFields.push_back(ptr);
  return 0;
}

// -----------------------------------------------------------------------------
// Tabulation utilities (slice / ranges / grids)
// -----------------------------------------------------------------------------
/**
 * @details With two sets (pos1, pos2), all unique layer‑pair combinations are
 * produced. If pos2 has a single point, three ghost points are added
 * (x,y,z offsets) to stabilize range estimation near singular pairs.
 */
std::vector<PositionGroup> DomainLayered3D::slice(const std::vector<rvec>& pos1,
                                                  std::vector<rvec>* pos2_opt) {
    std::vector<PositionGroup> result;

    if (pos2_opt == nullptr) {
        // Single set: group by layer index of pos1
        const auto& pos = pos1;
        std::map<int, PositionGroup> groups;
        for (size_t i = 0; i < pos.size(); ++i) {
            int layer_index = this->GetLayerIndex(pos[i]);
            if (groups.find(layer_index) == groups.end()) {
                groups[layer_index] = PositionGroup{static_cast<size_t>(layer_index), 0, {}, {}, {}, {}};
            }
            groups[layer_index].pos1.push_back(pos[i]);
            groups[layer_index].ind1.push_back(i);
        }
        for (std::map<int, PositionGroup>::iterator it = groups.begin(); it != groups.end(); ++it) {
            result.push_back(it->second);
        }
    } else {
      // Two sets of positions
      auto& pos2 = *pos2_opt;
      if (pos2.size() == 1) { // add three small displacements to avoid degeneracy
        pos2.push_back(pos2[0] + rvec(5e-4, 0.0, 0.0));
        pos2.push_back(pos2[0] + rvec(0.0, 5e-4, 0.0));
        pos2.push_back(pos2[0] + rvec(0.0, 0.0, 5e-4));
      }

      // Precompute unique layer indices for pos1 and pos2
      std::map<int, std::vector<size_t>> layer_map1, layer_map2;
      for (size_t i = 0; i < pos1.size(); ++i) {
          int layer_index = this->GetLayerIndex(pos1[i]);
          layer_map1[layer_index].push_back(i);
      }
      for (size_t j = 0; j < pos2.size(); ++j) {
          int layer_index = this->GetLayerIndex(pos2[j]);
          layer_map2[layer_index].push_back(j);
      }

      // Loop over unique layer indices and group positions
      for (const auto& layer1 : layer_map1) {
        int k1 = layer1.first;
        for (const auto& layer2 : layer_map2) {
          int k2 = layer2.first;
          const auto& indices1 = layer1.second;
          const auto& indices2 = layer2.second;

          PositionGroup group;
          group.i1 = k1;
          group.i2 = k2;

          // Group positions from pos1
          for (size_t i : indices1) {
              group.pos1.push_back(pos1[i]);
              group.ind1.push_back(i);
          }

          // Group positions from pos2
          for (size_t j : indices2) {
              group.pos2.push_back(pos2[j]);
              group.ind2.push_back(j);
          }
          result.push_back(group);
        }
      }
    }
    return result;
}

/**
 * @brief Adjusts the Z-coordinate range with a margin while ensuring boundaries.
 *
 * @details This helper function modifies the input Z-range by adding a margin while ensuring
 * that the adjusted range stays within the specified layer boundaries.
 *
 * @param z The original Z-coordinate range as a vector of two values (min, max).
 * @param zlayer The Z-coordinate boundaries of the layer.
 * @param zmin Minimum allowed vertical range.
 * @param margin Fraction of the range used to adjust the boundaries.
 * @return A vector containing the adjusted Z-range (min, max).
 */
std::vector<double> adjustZRange(const std::vector<double>& z,
                                 const std::vector<double>& zlayer,
                                 double zmin, double margin) {
    double dz = margin * (z[1] - z[0]);
    return {
        std::max(z[0] - dz, zlayer[0] + zmin),
        std::min(z[1] + dz, zlayer[1] - zmin)
    };
}

/**
 * @details This function determines the spatial range in both the radial and vertical (Z)
 * directions for a given set of points in a layered medium. It ensures that the
 * computed range remains within valid boundaries while applying necessary adjustments
 * using a margin parameter.
 */
Range DomainLayered3D::calculateRange(const PositionGroup& pts, double rmin,
                                      double zmin, double margin) {
  Range range;
  range.i1 = pts.i1;
  range.i2 = pts.i2;

  // Compute radial distance range
  double dmin = std::numeric_limits<double>::max();
  double dmax = std::numeric_limits<double>::lowest();
  for (const auto& p1 : pts.pos1) {
    for (const auto& p2 : pts.pos2) {
      double dx = p1[0] - p2[0];
      double dy = p1[1] - p2[1];
      double distance = std::sqrt(dx * dx + dy * dy);
      dmin = std::min(dmin, distance);
      dmax = std::max(dmax, distance);
    }
  }
  if (pts.ind2.size() == 1) {
    // Points in the same layer
    range.r = {dmin - rmin, dmin + rmin};
  } else {
    // Points in different layers
    range.r = {std::max(dmin, rmin), (1.0 + margin) * dmax};
  }

  // Compute z-ranges
  std::vector<double> z1 = {std::numeric_limits<double>::max(),
                            std::numeric_limits<double>::lowest()};
  std::vector<double> z2 = z1;

  for (const auto& p1 : pts.pos1) {
    z1[0] = std::min(z1[0], p1[2]);
    z1[1] = std::max(z1[1], p1[2]);
  }

  for (const auto& p2 : pts.pos2) {
    z2[0] = std::min(z2[0], p2[2]);
    z2[1] = std::max(z2[1], p2[2]);
  }

  if (z1[0] == z1[1]) z1[1] += zmin;
  if (z2[0] == z2[1]) z2[1] += zmin;

  if (range.i1 == range.i2) {
    // Points in the same layer
    if (range.i1 == 0) {
      // Lowest layer
      std::vector<double> zlim = {
        2 * this->zValsInterfaces[0] - (z1[1] + z2[1]),
        2 * this->zValsInterfaces[0] - (z1[0] + z2[0])
      };
      range.z1 = adjustZRange(zlim, {0, 1e30}, zmin, 2 * margin);
      range.z2.clear();
    } else if (range.i1 == this->zValsInterfaces.size()) {
      // Highest layer
      std::vector<double> zlim = {
        (z1[0] + z2[0]) - 2 * this->zValsInterfaces.back(),
        (z1[1] + z2[1]) - 2 * this->zValsInterfaces.back()
      };
      range.z1 = adjustZRange(zlim, {0, 1e30}, zmin, 2 * margin);
      range.z2.clear();
    } else {
      // Inside the layer structure
      z1 = adjustZRange(z1, {this->zValsInterfaces[range.i1 - 1], this->zValsInterfaces[range.i1]}, zmin, 2 * margin);
      z2 = adjustZRange(z2, {this->zValsInterfaces[range.i1 - 1], this->zValsInterfaces[range.i1]}, zmin, 2 * margin);

      range.z1 = {
        z1[0] - z2[1] + (this->zValsInterfaces[range.i1] - this->zValsInterfaces[range.i1 - 1]),
        z1[1] - z2[0] + (this->zValsInterfaces[range.i1] - this->zValsInterfaces[range.i1 - 1])
      };
      range.z2 = {
        z1[0] + z2[0] - 2 * this->zValsInterfaces[range.i1 - 1],
        z1[1] + z2[1] - 2 * this->zValsInterfaces[range.i1 - 1]
      };
    }
  }
  else {
    // Points in different layers
    std::vector<double> zTemp = {-1e30};
    zTemp.insert(zTemp.end(), this->zValsInterfaces.begin(), this->zValsInterfaces.end());
    zTemp.push_back(1e30);
    range.z1 = adjustZRange(z1, {zTemp[range.i1], zTemp[range.i1 + 1]}, zmin, 2 * margin);
    range.z2 = adjustZRange(z2, {zTemp[range.i2], zTemp[range.i2 + 1]}, zmin, 2 * margin);
  }
  return range;
}

/**
 * @details This function generates a structured grid in the radial and vertical (Z) directions
 * for each position group, ensuring appropriate spacing and range adjustments. It uses
 * the calculated spatial bounds to define the grid points.
 */
std::vector<Grid> DomainLayered3D::computeGrid(const std::vector<PositionGroup>& pos,
                                               int nr, int nz, double rmin,
                                               double margin, double zmin) {
    std::vector<Grid> grids;

    // Compute ranges for all position groups
    for (const auto& group : pos) {
        Range range = calculateRange(group, rmin, zmin, margin);
        Grid grid;
        grid.i1 = range.i1;
        grid.i2 = range.i2;
        int nr_ = nr, nz_ = nz;
        if (!range.z2.empty()) { nr_ = nr / 2; nz_ = nz / 2; }

        // Compute radial grid
        double r_start = range.r[0];
        double r_end = range.r[1];
        grid.r.resize(nr_);
        for (int i = 0; i < nr_; ++i) {
            grid.r[i] = r_start + i * (r_end - r_start) / (nr_ - 1);
        }

        // Compute Z1 grid
        double z1_start = range.z1[0];
        double z1_end = range.z1[1];
        grid.z1.resize(nz_);
        for (int i = 0; i < nz_; ++i) {
            grid.z1[i] = z1_start + i * (z1_end - z1_start) / (nz_ - 1);
        }

        // Compute Z2 grid if it exists
        if (!range.z2.empty()) {
            double z2_start = range.z2[0];
            double z2_end = range.z2[1];
            grid.z2.resize(nz_);
            for (int i = 0; i < nz_; ++i) {
                grid.z2[i] = z2_start + i * (z2_end - z2_start) / (nz_ - 1);
            }
        } else {
            grid.z2.clear();
        }

        grids.push_back(grid);
    }

    return grids;
}

// -----------------------------------------------------------------------------
// Fields and cross sections (layered post‑processing)
// -----------------------------------------------------------------------------
cvec DomainLayered3D::SecFieldE(blitz::Array<dcmplx, 1> *solVector, rvec pos) {
  cvec result(0, 0, 0);
  int offset(solVector->shape()(0) / 2);
  using std::vector;
  using std::pair;
  vector<RWGFun>::iterator rwgM, rwgN;

  // Partition the RWG functions to groups over same triangle
  vector<RWGFun *> rwgPerTri;
  for (rwgM = rwgFuns.begin(); rwgM != rwgFuns.end(); rwgM = rwgN) {
    rwgPerTri.clear();
    for (rwgN = rwgM;
         rwgN != rwgFuns.end() and rwgN->TrianglePtr() == rwgM->TrianglePtr();
         rwgN++) {
      rwgPerTri.push_back(&(*rwgN));
    }
    vector<pair<cvec, cvec> > DeltaKappaHom(
        rwgPerTri.size(), make_pair(cvec(0., 0., 0.), cvec(0., 0., 0.)));

     vector<pair<cvec, cvec> > DeltaKappaE(
        rwgPerTri.size(), make_pair(cvec(0., 0., 0.), cvec(0., 0., 0.)));

    rvec jCenter = rwgPerTri[0]->TrianglePtr()->Center();
    if (GetLayerIndex(pos) == GetLayerIndex(jCenter)) {
      grnFun[GetLayerIndex(pos)].SameTriDeltaKappa(pos, rwgPerTri, DeltaKappaHom);
    }
    grnFunLayeredFields->SameTriDeltaKappa(pos, rwgPerTri, DeltaKappaE);
    
    for (int i = 0; i < (int)rwgPerTri.size(); ++i) {
      result -= I * omega * MU0 * Mu(jCenter) * 
                (DeltaKappaHom[i].first + DeltaKappaE[i].first) *
                (*solVector)(rwgPerTri[i]->EdgeIndex()) +
                (DeltaKappaHom[i].second + DeltaKappaE[i].second) *
                (*solVector)(rwgPerTri[i]->EdgeIndex() + offset);
    }
  }
  return result;
}

cvec DomainLayered3D::IncFieldE(rvec pos) {
  cvec incField(0, 0, 0);
  for (std::vector<IncidentField *>::iterator incIter = incidentFields.begin();
       incIter != incidentFields.end(); incIter++) {
    std::vector<cvec> incFieldVec = (*incIter)->EvaluateLayeredEandH(pos);
    incField += incFieldVec[0]; // The first element is E field
  }
  return incField;
}

cvec DomainLayered3D::SecFieldH(blitz::Array<dcmplx, 1> *solVector, rvec pos) {
  cvec result(0, 0, 0);
  int offset(solVector->shape()(0) / 2);
  using std::vector;
  using std::pair;
  vector<RWGFun>::iterator rwgM, rwgN;

  // Partition the RWG functions to groups over same triangle
  vector<RWGFun *> rwgPerTri;
  for (rwgM = rwgFuns.begin(); rwgM != rwgFuns.end(); rwgM = rwgN) {
    rwgPerTri.clear();
    for (rwgN = rwgM;
         rwgN != rwgFuns.end() and rwgN->TrianglePtr() == rwgM->TrianglePtr();
         rwgN++) {
      rwgPerTri.push_back(&(*rwgN));
    }
    vector<pair<cvec, cvec> > DeltaKappaHom(
        rwgPerTri.size(), make_pair(cvec(0., 0., 0.), cvec(0., 0., 0.)));

     vector<pair<cvec, cvec> > DeltaKappaH(
        rwgPerTri.size(), make_pair(cvec(1., 0., 0.), cvec(0., 0., 0.)));

    rvec jCenter = rwgPerTri[0]->TrianglePtr()->Center();
    if (GetLayerIndex(pos) == GetLayerIndex(jCenter)) {
      grnFun[GetLayerIndex(pos)].SameTriDeltaKappa(pos, rwgPerTri, DeltaKappaHom);
    }
    grnFunLayeredFields->SameTriDeltaKappa(pos, rwgPerTri, DeltaKappaH);
    
    for (int i = 0; i < (int)rwgPerTri.size(); ++i) {
      result -= I * omega * EPS0 * Epsilon(jCenter) * 
                (DeltaKappaHom[i].first + DeltaKappaH[i].first) *
                (*solVector)(rwgPerTri[i]->EdgeIndex() + offset) -
                (DeltaKappaHom[i].second + DeltaKappaH[i].second) * 
                (*solVector)(rwgPerTri[i]->EdgeIndex());
    }
  }
  return result;
}

std::pair<cvec, cvec> DomainLayered3D::SecFieldEandH(blitz::Array<dcmplx, 1> *solVector,
                                                     rvec pos) {
  std::pair<cvec, cvec> result(make_pair(cvec(0, 0, 0), cvec(0, 0, 0)));
  int offset(solVector->shape()(0) / 2);
  using std::vector;
  using std::pair;
  vector<RWGFun>::iterator rwgM, rwgN;

  // Partition the RWG functions to groups over same triangle
  vector<RWGFun *> rwgPerTri;
  for (rwgM = rwgFuns.begin(); rwgM != rwgFuns.end(); rwgM = rwgN) {
    rwgPerTri.clear();
    for (rwgN = rwgM;
         rwgN != rwgFuns.end() and rwgN->TrianglePtr() == rwgM->TrianglePtr();
         rwgN++) {
      rwgPerTri.push_back(&(*rwgN));
    }
    vector<pair<cvec, cvec> > DeltaKappaHom(
        rwgPerTri.size(), make_pair(cvec(0., 0., 0.), cvec(0., 0., 0.)));

     vector<pair<cvec, cvec> > DeltaKappaE(
        rwgPerTri.size(), make_pair(cvec(0., 0., 0.), cvec(0., 0., 0.)));

     vector<pair<cvec, cvec> > DeltaKappaH(
        rwgPerTri.size(), make_pair(cvec(1., 0., 0.), cvec(0., 0., 0.)));
    
    rvec jCenter = rwgPerTri[0]->TrianglePtr()->Center();
    if (GetLayerIndex(pos) == GetLayerIndex(jCenter)) {
      grnFun[GetLayerIndex(pos)].SameTriDeltaKappa(pos, rwgPerTri, DeltaKappaHom);
    }
    grnFunLayeredFields->SameTriDeltaKappa(pos, rwgPerTri, DeltaKappaE);
    grnFunLayeredFields->SameTriDeltaKappa(pos, rwgPerTri, DeltaKappaH);

    for (int i = 0; i < (int)rwgPerTri.size(); ++i) {
      result.first -= I * omega * MU0 * Mu(jCenter) * 
                      (DeltaKappaHom[i].first + DeltaKappaE[i].first) *
                      (*solVector)(rwgPerTri[i]->EdgeIndex()) +
                      (DeltaKappaHom[i].second + DeltaKappaE[i].second) *
                      (*solVector)(rwgPerTri[i]->EdgeIndex() + offset);
      result.second -= I * omega * EPS0 * Epsilon(jCenter) * 
                       (DeltaKappaHom[i].first + DeltaKappaH[i].first) *
                       (*solVector)(rwgPerTri[i]->EdgeIndex() + offset) -
                       (DeltaKappaHom[i].second + DeltaKappaH[i].second) * 
                       (*solVector)(rwgPerTri[i]->EdgeIndex());
    }
  }
#ifdef DEBUG
  dumpField("Efield", result.first);
  dumpField("Hfield", result.second);
#endif
  return result;
}

cvec DomainLayered3D::IncFieldH(rvec pos) {
  cvec incField(0, 0, 0);
  for (std::vector<IncidentField *>::iterator incIter = incidentFields.begin();
       incIter != incidentFields.end(); incIter++) {
    std::vector<cvec> incFieldVec = (*incIter)->EvaluateLayeredEandH(pos);
    incField += incFieldVec[1]; // The second element is H field
  }
  return incField;
}

/**
 * @brief Compute the element-wise complex conjugate of a 3-vector.
 * @param in Input complex vector.
 * @return Complex vector containing conjugated components of @p in.
 */
cvec conjVecL(cvec in) { return cvec(conj(in(0)), conj(in(1)), conj(in(2))); }

/**
 * @brief Extract the real part of each component of a complex 3-vector.
 * @param in Input complex vector.
 * @return Real vector with the real parts of @p in.
 */
rvec realVecL(cvec in) { return rvec(real(in(0)), real(in(1)), real(in(2))); }

std::vector<double> DomainLayered3D::scatteringCS(blitz::Array<dcmplx,1>* solVector)
{
  // Gather unique triangles & normals
  std::vector<Triangle*> vTri;
  std::vector<rvec> vNormals;
  vTri.reserve(rwgFuns.size());
  vNormals.reserve(rwgFuns.size());

  for (RWGFun& rwg : rwgFuns) {
    Triangle* T = rwg.TrianglePtr();
    auto it = std::find(vTri.begin(), vTri.end(), T);
    if (it == vTri.end()) {
      vTri.push_back(T);
      vNormals.push_back(T->FBNormalVector(rwg.Side()));
    }
  }

  size_t Ntri = vTri.size();  // number of unique triangles

  // Quadrature rule on reference triangle
  int nQuad = dunavant::N5;    // number of quadrature points
  auto& xQuad = dunavant::x5;  // barycentrics [nQuad][3]
  auto& wQuad = dunavant::w5;  // weights[nQuad]

  // Allocate per-triangle, per-quad-point accumulators
  std::vector<std::vector<cvec>> Etot(Ntri, std::vector<cvec>(nQuad, cvec(0, 0, 0))),
                                 Htot(Ntri, std::vector<cvec>(nQuad, cvec(0, 0, 0)));

  int offset = solVector->shape()(0) / 2; // offset for β coefficients

  // Fill E, H normal and tangential components at each sample point
  for (RWGFun& rwg : rwgFuns) {
    // find which triangle index
    Triangle* T = rwg.TrianglePtr();
    size_t triIndex = std::distance(vTri.begin(),
                      std::find(vTri.begin(), vTri.end(), T));
    rvec normal = vNormals[triIndex];

    // material parameters (step off interface if needed)
    rvec ctr = T->Center();
    bool onInt = false;
    for (double z0 : zValsInterfaces)
      if (std::abs( ctr[2] - z0 ) < 1e-6) { onInt = true; break; }

    // element local material parameters
    dcmplx epsLoc = onInt ? Epsilon(ctr + 1e-1 * normal) : Epsilon(ctr);
    dcmplx muLoc  = onInt ? Mu(ctr + 1e-1 * normal) : Mu(ctr);

    // RWG coefficients
    dcmplx alpha = (*solVector)(rwg.EdgeIndex());
    dcmplx beta  = (*solVector)(rwg.EdgeIndex() + offset);
    double divf  = rwg.RWGPreFactor();

    // triangle nodes
    rvec r[3] = { T->Node(0), T->Node(1), T->Node(2) };

    // loop quad points
    for (int q = 0; q < nQuad; ++q) {
      // physical sample point
      rvec s = xQuad[q][0] * r[0] + xQuad[q][1] * r[1] + xQuad[q][2] * r[2];

      // shape‐function at s
      rvec fval = rwg.Evaluate(s);

      // normal‐component (constant over tri)
      cvec En_contrib = -I * alpha * divf * normal / (omega * epsLoc * EPS0);
      cvec Hn_contrib = -I * beta * divf * normal / (omega * muLoc * MU0);

      // tangential‐component
      cvec Et_contrib =  cross( (cvec)normal, (cvec)(beta * fval));
      cvec Ht_contrib = -cross( (cvec)normal, (cvec)(alpha * fval));

      // J,M currents have contributions from each RWG with opposite sign, thus -=
      Etot[triIndex][q] -= En_contrib + Et_contrib;
      Htot[triIndex][q] -= Hn_contrib + Ht_contrib;
    }
  }

  // Quadrature - find the scattering cross section
  double sigma = 0.0;
  double near = 0.0;
  for (size_t t = 0; t < Ntri; ++t) {
    rvec normal = vNormals[t];
    double area = vTri[t]->Area();
    
    // re-compute nodes for this triangle
    rvec r[3] = { vTri[t]->Node(0), vTri[t]->Node(1), vTri[t]->Node(2) };

    for (int q = 0; q < nQuad; ++q) {
      // map barycentrics -> point
      rvec s = xQuad[q][0] * r[0] + xQuad[q][1] * r[1] + xQuad[q][2] * r[2];

      cvec Einc = IncFieldE(s);
      cvec Hinc = IncFieldH(s);
      rvec Sinc = realVecL( cross(Einc, conjVecL(Hinc)) ) / 2.0;

      cvec Escat = Etot[t][q] - Einc;
      cvec Hscat = Htot[t][q] - Hinc;
      rvec Sscat = realVecL( cross(Escat, conjVecL(Hscat)) ) / 2.0;

      sigma += wQuad[q] * area * dot(normal, Sscat) / sqrt( dot( Sinc, Sinc ) );
      near += real( wQuad[q] * area * dot( Escat, conjVecL(Escat) ) );
    }
  }
  return { sigma, near };
}

double DomainLayered3D::absorptionCS(blitz::Array<dcmplx,1>* solVector)
{
  // Gather unique triangles & normals
  std::vector<Triangle*> vTri;
  std::vector<rvec> vNormals;
  vTri.reserve(rwgFuns.size());
  vNormals.reserve(rwgFuns.size());

  for (RWGFun& rwg : rwgFuns) {
    Triangle* T = rwg.TrianglePtr();
    auto it = std::find(vTri.begin(), vTri.end(), T);
    if (it == vTri.end()) {
      vTri.push_back(T);
      vNormals.push_back(T->FBNormalVector(rwg.Side()));
    }
  }

  size_t Ntri = vTri.size();  // number of unique triangles

  // Quadrature rule on reference triangle
  int nQuad = dunavant::N5;    // number of quadrature points
  auto& xQuad = dunavant::x5;  // barycentrics [nQuad][3]
  auto& wQuad = dunavant::w5;  // weights[nQuad]

  // Allocate per-triangle, per-quad-point accumulators
  std::vector<std::vector<cvec>> Etot(Ntri, std::vector<cvec>(nQuad, cvec(0, 0, 0))),
                                 Htot(Ntri, std::vector<cvec>(nQuad, cvec(0, 0, 0)));

  int offset = solVector->shape()(0) / 2; // offset for β coefficients

  // Fill E, H normal and tangential components at each sample point
  for (RWGFun& rwg : rwgFuns) {
    // find which triangle index
    Triangle* T = rwg.TrianglePtr();
    size_t triIndex = std::distance(vTri.begin(),
                      std::find(vTri.begin(), vTri.end(), T));
    rvec normal = vNormals[triIndex];

    // material parameters (step off interface if needed)
    rvec ctr = T->Center();
    bool onInt = false;
    for (double z0 : zValsInterfaces)
      if (std::abs( ctr[2] - z0 ) < 1e-6) { onInt = true; break; }

    // element local material parameters
    dcmplx epsLoc = onInt ? Epsilon(ctr + 1e-1 * normal) : Epsilon(ctr);
    dcmplx muLoc  = onInt ? Mu(ctr + 1e-1 * normal) : Mu(ctr);

    // RWG coefficients
    dcmplx alpha = (*solVector)(rwg.EdgeIndex());
    dcmplx beta  = (*solVector)(rwg.EdgeIndex() + offset);
    double divf  = rwg.RWGPreFactor();

    // triangle nodes
    rvec r[3] = { T->Node(0), T->Node(1), T->Node(2) };

    // loop quad points
    for (int q = 0; q < nQuad; ++q) {
      // physical sample point
      rvec s = xQuad[q][0] * r[0] + xQuad[q][1] * r[1] + xQuad[q][2] * r[2];

      // shape‐function at s
      rvec fval = rwg.Evaluate(s);

      // normal‐component (constant over tri)
      cvec En_contrib = -I * alpha * divf * normal / (omega * epsLoc * EPS0);
      cvec Hn_contrib = -I * beta * divf * normal / (omega * muLoc * MU0);

      // tangential‐component
      cvec Et_contrib =  cross( (cvec)normal, (cvec)(beta * fval));
      cvec Ht_contrib = -cross( (cvec)normal, (cvec)(alpha * fval));

      // J,M currents have contributions from each RWG with opposite sign, thus -=
      Etot[triIndex][q] -= En_contrib + Et_contrib;
      Htot[triIndex][q] -= Hn_contrib + Ht_contrib;
    }
  }

  // Quadrature - find the scattering cross section
  double sigma = 0.0;
  for (size_t t = 0; t < Ntri; ++t) {
    rvec normal = vNormals[t];
    double area = vTri[t]->Area();
    
    // re-compute nodes for this triangle
    rvec r[3] = { vTri[t]->Node(0), vTri[t]->Node(1), vTri[t]->Node(2) };

    for (int q = 0; q < nQuad; ++q) {
      // map barycentrics -> point
      rvec s = xQuad[q][0] * r[0] + xQuad[q][1] * r[1] + xQuad[q][2] * r[2];

      cvec Einc = IncFieldE(s);
      cvec Hinc = IncFieldH(s);
      rvec Sinc = realVecL( cross(Einc, conjVecL(Hinc)) ) / 2.0;

      cvec E = Etot[t][q];
      cvec H = Htot[t][q];
      rvec S = realVecL( cross(E, conjVecL(H)) ) / 2.0;

      sigma -= wQuad[q] * area * dot(normal, S) / sqrt( dot( Sinc, Sinc ) );
    }
  }
  return sigma;
}

double DomainLayered3D::extinctionCS(blitz::Array<dcmplx,1>* solVector)
{
  // Gather unique triangles & normals
  std::vector<Triangle*> vTri;
  std::vector<rvec> vNormals;
  vTri.reserve(rwgFuns.size());
  vNormals.reserve(rwgFuns.size());

  for (RWGFun& rwg : rwgFuns) {
    Triangle* T = rwg.TrianglePtr();
    auto it = std::find(vTri.begin(), vTri.end(), T);
    if (it == vTri.end()) {
      vTri.push_back(T);
      vNormals.push_back(T->FBNormalVector(rwg.Side()));
    }
  }

  size_t Ntri = vTri.size();  // number of unique triangles

  // Quadrature rule on reference triangle
  int nQuad = dunavant::N5;    // number of quadrature points
  auto& xQuad = dunavant::x5;  // barycentrics [nQuad][3]
  auto& wQuad = dunavant::w5;  // weights[nQuad]

  // Allocate per-triangle, per-quad-point accumulators
  std::vector<std::vector<cvec>> Etot(Ntri, std::vector<cvec>(nQuad, cvec(0, 0, 0))),
                                 Htot(Ntri, std::vector<cvec>(nQuad, cvec(0, 0, 0)));

  int offset = solVector->shape()(0) / 2; // offset for β coefficients

  // Fill E, H normal and tangential components at each sample point
  for (RWGFun& rwg : rwgFuns) {
    // find which triangle index
    Triangle* T = rwg.TrianglePtr();
    size_t triIndex = std::distance(vTri.begin(),
                      std::find(vTri.begin(), vTri.end(), T));
    rvec normal = vNormals[triIndex];

    // material parameters (step off interface if needed)
    rvec ctr = T->Center();
    bool onInt = false;
    for (double z0 : zValsInterfaces)
      if (std::abs( ctr[2] - z0 ) < 1e-6) { onInt = true; break; }

    // element local material parameters
    dcmplx epsLoc = onInt ? Epsilon(ctr + 1e-1 * normal) : Epsilon(ctr);
    dcmplx muLoc  = onInt ? Mu(ctr + 1e-1 * normal) : Mu(ctr);

    // RWG coefficients
    dcmplx alpha = (*solVector)(rwg.EdgeIndex());
    dcmplx beta  = (*solVector)(rwg.EdgeIndex() + offset);
    double divf  = rwg.RWGPreFactor();

    // triangle nodes
    rvec r[3] = { T->Node(0), T->Node(1), T->Node(2) };

    // loop quad points
    for (int q = 0; q < nQuad; ++q) {
      // physical sample point
      rvec s = xQuad[q][0] * r[0] + xQuad[q][1] * r[1] + xQuad[q][2] * r[2];

      // shape‐function at s
      rvec fval = rwg.Evaluate(s);

      // normal‐component (constant over tri)
      cvec En_contrib = -I * alpha * divf * normal / (omega * epsLoc * EPS0);
      cvec Hn_contrib = -I * beta * divf * normal / (omega * muLoc * MU0);

      // tangential‐component
      cvec Et_contrib =  cross( (cvec)normal, (cvec)(beta * fval));
      cvec Ht_contrib = -cross( (cvec)normal, (cvec)(alpha * fval));

      // J,M currents have contributions from each RWG with opposite sign, thus -=
      Etot[triIndex][q] -= En_contrib + Et_contrib;
      Htot[triIndex][q] -= Hn_contrib + Ht_contrib;
    }
  }

  // Quadrature - find the scattering cross section
  double sigma = 0.0;
  for (size_t t = 0; t < Ntri; ++t) {
    rvec normal = vNormals[t];
    double area = vTri[t]->Area();
    
    // re-compute nodes for this triangle
    rvec r[3] = { vTri[t]->Node(0), vTri[t]->Node(1), vTri[t]->Node(2) };

    for (int q = 0; q < nQuad; ++q) {
      // map barycentrics -> point
      rvec s = xQuad[q][0] * r[0] + xQuad[q][1] * r[1] + xQuad[q][2] * r[2];

      cvec Einc = IncFieldE(s);
      cvec Hinc = IncFieldH(s);
      rvec Sinc = realVecL( cross(Einc, conjVecL(Hinc)) ) / 2.0;

      cvec Escat = Etot[t][q] - Einc;
      cvec Hscat = Htot[t][q] - Hinc;
      rvec Sext  = realVecL( cross(Einc, conjVecL(Hscat)) + 
                             cross(Escat, conjVecL(Hinc)) ) / 2.0;

      sigma -= wQuad[q] * area * dot(normal, Sext) / sqrt( dot( Sinc, Sinc ) );
    }
  }
  return sigma;
}

// -----------------------------------------------------------------------------
// Debugging and testing
// -----------------------------------------------------------------------------
void DomainLayered3D::test_slice()
{
  std::string filename = "positions.txt";
  std::vector<rvec> points = readPositionsFromFile(filename);

  // Example usage: Single set of positions
  auto groups1 = slice(points);
  std::cout << std::fixed << std::setprecision(2); // Set precision to 2 digits
  std::cout << "\n" << "----------------------------------------------------------------------\n";
  std::cout << "\n" << "            This is the tabulation slice calculation test:\n";
  std::cout << "\n" << "----------------------------------------------------------------------\n";
  std::cout << "Single set of positions:\n";
  for (const auto& group : groups1) {
      std::cout << "Layer " << group.i1 << ": ";
      for (const auto& pos : group.pos1) {
          std::cout << "(" << pos[0] << ", " << pos[1] << ", " << pos[2] << "),\n";
      }
      std::cout << "\n";
  }

  // Example usage: Two sets of positions
  auto groups2 = slice(points, &points);
  std::cout << "\nTwo sets of positions:\n";
  for (const auto& group : groups2) {
    std::cout << "Layer pair (" << group.i1 << ", " << group.i2 << "):\n";
    std::cout << "--Positions from set 1: ";
    std::cout << group.pos1.size() << "\n";
    std::cout << "--Indices from set 1: ";
    std::cout << group.ind1.size() << "\n";
    std::cout << "\n" << "--Positions from set 2: ";
    std::cout << group.pos2.size() << "\n";
    std::cout << "--Indices from set 2: ";
    std::cout << group.ind2.size() << "\n";
    std::cout << "\n\n";
  }
}

void DomainLayered3D::test_range() 
{
  std::string filename = "positions.txt";
  std::vector<rvec> points = readPositionsFromFile(filename);

  auto groups1 = slice(points, &points);
  std::vector<Range> ranges;
  std::cout << std::fixed << std::setprecision(6); // Set precision to 6 digits
  std::cout << "\n" << "----------------------------------------------------------------------\n";
  std::cout << "\n" << "            This is the tabulation range calculation test:\n";
  std::cout << "\n" << "----------------------------------------------------------------------\n";
  for (const auto& pt : groups1) {
    Range range = calculateRange(pt);
    ranges.push_back(range);
    std::cout << "Layer pair (" << range.i1 << ", " << range.i2 << "): ";
    std::cout << "Radial range: [" << range.r[0] << ", " << range.r[1] << "], ";
    std::cout << "Z1 range: [" << range.z1[0] << ", " << range.z1[1] << "], ";
    std::cout << "Z2 range: [" << range.z2[0] << ", " << range.z2[1] << "]\n";
  }
  std::cout << "\n";
}

void DomainLayered3D::test_grid() 
{
  std::string filename = "positions.txt";
  std::vector<rvec> points = readPositionsFromFile(filename);

  auto groups1 = slice(points, &points);
  std::vector<Grid> grids = computeGrid(groups1, 5, 5);
  std::cout << std::fixed << std::setprecision(5); // Set precision to 5 digits
  std::cout << "\n" << "----------------------------------------------------------------------\n";
  std::cout << "\n" << "           This is the tabulation grid calculation test:\n";
  std::cout << "\n" << "----------------------------------------------------------------------\n";
  // Print grids
  for (size_t i = 0; i < grids.size(); ++i) {
    std::cout << "Grid " << i + 1 << ":\n";
    std::cout << "  r grid: ";
    for (const auto& r : grids[i].r) {
      std::cout << r << " ";
    }
    std::cout << "\n  Z1 grid: ";
    for (const auto& z1 : grids[i].z1) {
      std::cout << z1 << " ";
    }
    if (!grids[i].z2.empty()) {
      std::cout << "\n  Z2 grid: ";
      for (const auto& z2 : grids[i].z2) {
        std::cout << z2 << " ";
      }
    }
    std::cout << "\n\n";
  }
  std::cout << "\n";
}

void DomainLayered3D::dumpField(const std::string& fieldName,
                                const blitz::TinyVector<std::complex<double>, 3>& vec) {
    static std::ofstream eFile("Efield.txt");
    static std::ofstream hFile("Hfield.txt");

    std::ofstream* out = nullptr;
    if (fieldName == "Efield") {
        out = &eFile;
    } else if (fieldName == "Hfield") {
        out = &hFile;
    } else {
        throw std::runtime_error("Unsupported field name: " + fieldName);
    }

    *out << std::fixed << std::setprecision(6);
    for (int i = 0; i < 3; ++i) {
        const auto& c = vec(i);
        *out << c.real();
        if (c.imag() >= 0) *out << " + ";
        else *out << " - ";
        *out << std::abs(c.imag()) << "i";
        if (i < 2) *out << " ";
    }
    *out << "; ";
}

/**
 * @details This method first calls the base class implementation to clear the 
 * @f$ \text{incidentFields} @f$ list (deleting the Dipole objects). It then 
 * deletes and nullifies the @f$ \text{grnFunLayeredFields} @f$ pointer.
 * This is crucial because Dipole objects are recreated for each subjob, 
 * but the Domain object persists. Resetting @f$ \text{grnFunLayeredFields} @f$ 
 * forces @f$ \text{SimJob} @f$ to call @f$ \text{newGrnFunLayered()} @f$ 
 * again, which correctly wires the new Dipole sources.
 */
int DomainLayered3D::ClearIncidentFields() {
  Domain::ClearIncidentFields(); 
  if (grnFunLayeredFields) {
      delete grnFunLayeredFields;
      grnFunLayeredFields = nullptr;
  }
  return 0;
}
