/**
 * @file LayeredMediaUtils.cpp
 * @brief Implementation for the LayeredMediaUtils class.
 */

#include "LayeredMediaUtils.h"
#include <algorithm>
#include <numeric>
#include <unordered_map>
#include <unordered_set>

/**
 * @brief Computes the square root of a complex number with proper handling of negative imaginary parts.
 *
 * This function ensures correct computation of the square root for cases where the imaginary part 
 * is negative or zero to maintain consistency in branch cuts.
 *
 * @param arg The input complex number.
 * @return The square root of the input complex number.
 */
dcmplx complexSqrt(dcmplx arg) {
  if (imag(arg) == 0) {
    dcmplx buf(real(arg), 0);
    arg = buf;
  }
  if (imag(arg) < 0) {
    return -sqrt(arg);
  } else {
    return sqrt(arg);
  }
}

/**
 * @details This constructor initializes the wavenumbers, permittivity, permeability, 
 * interface positions, and layer thicknesses for a layered medium.
 */
LayeredMediaUtils::LayeredMediaUtils(const std::vector<dcmplx>& inLayeredK,
                                     const std::vector<dcmplx>& inEpsilonMedium,
                                     const std::vector<dcmplx>& inMuMedium,
                                     const std::vector<double>& inZValsInterfaces,
                                     const std::vector<double>& inThickness,
                                     int inMidLayerIndex) 
    : k_L(inLayeredK),
      epsilon(inEpsilonMedium),
      mu(inMuMedium),
      zValsInterfaces(inZValsInterfaces),
      thickness(inThickness),
      midLayerIndex(inMidLayerIndex) {}

/**
 * @details This function identifies the index of the layer in which a given point resides 
 * based on its z-coordinate relative to the predefined interface positions.
 */
int LayeredMediaUtils::GetLayerIndex(rvec pos) const {
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

/**
 * @details This function retrieves the Fresnel coefficients for the interface between two given 
 * layers, considering the specified transverse wavenumber and wave polarization (TE or TM).
 */
std::vector<dcmplx> LayeredMediaUtils::FresnelCoeff(int layerIndexI, int layerIndexJ, 
                                                    dcmplx krho, const std::string& waves) const {
  
  // Validate layer indices
  if (layerIndexI < 0 || layerIndexI >= epsilon.size() || 
      layerIndexJ < 0 || layerIndexJ >= epsilon.size() || 
      std::abs(layerIndexI - layerIndexJ) > 1){
    throw std::out_of_range("Invalid layer indeces for Fresnel reflection coefficient calculation.");
  }

  // Check if the wave type is valid (TE or TM)
  if (waves != "TE" && waves != "TM"){
    throw std::invalid_argument("Invalid wave type. Use 'TE' or 'TM'.");
  }

  // Get the permittivity (TM) and permeability (TE) of the two layers
  dcmplx pI(0,0);
  dcmplx pJ(0,0);
  if (waves == "TM"){
    pI = epsilon[layerIndexI];
    pJ = epsilon[layerIndexJ];
  }
  else if (waves == "TE"){
    pI = mu[layerIndexI];
    pJ = mu[layerIndexJ];
  }
  
  // z component of the wavevector
  dcmplx kIz = k_L[layerIndexI];
  dcmplx kJz = k_L[layerIndexJ];
  kIz = kIz * kIz - krho * krho; kIz = complexSqrt(kIz);
  kJz = kJz * kJz - krho * krho; kJz = complexSqrt(kJz);

  // Fresnel reflection coefficient
  dcmplx fresnelReflCoeff;
  fresnelReflCoeff = (pJ * kIz - pI * kJz) /
                     (pJ * kIz + pI * kJz);

  // Fresnel transmission coefficient
  dcmplx fresnelTransCoeff;
  fresnelTransCoeff = (2.0 * pJ * kIz) / ( pJ * kIz + pI * kJz );

  std::vector<dcmplx> fresnelCoeffs = {fresnelReflCoeff, fresnelTransCoeff};
  
  return fresnelCoeffs;
}

/**
 * @details This function retrieves the transfer matrix that describes wave propagation 
 * across the interface between two specified layers, considering the given 
 * transverse wavenumber and wave polarization (TE or TM).
 */
std::vector<std::vector<dcmplx>> LayeredMediaUtils::InterfaceTransferMatrix(
    int layerIndexI, int layerIndexJ, dcmplx krho, const std::string& waves) const {

    // Validate layer indices
    if (layerIndexI < 0 || layerIndexI >= epsilon.size() || 
        layerIndexJ < 0 || layerIndexJ >= epsilon.size() || 
        std::abs(layerIndexI - layerIndexJ) > 1) {
        throw std::out_of_range("Invalid layer indices for transfer matrix calculation.");
    }

    // Validate wave type
    if (waves != "TE" && waves != "TM") {
        throw std::invalid_argument("Invalid wave type. Use 'TE' or 'TM'.");
    }

    // Calculate Fresnel reflection and transmission coefficients
    std::vector<dcmplx> fresnelCoeffs = this->FresnelCoeff(layerIndexI, layerIndexJ, krho, waves);
    dcmplx r = fresnelCoeffs[0];  // Reflection coefficient
    dcmplx t = fresnelCoeffs[1];  // Transmission coefficient

    // Prepare the transfer matrix for input polarization (TE or TM)
    std::vector<std::vector<dcmplx>> tMatrix(2, std::vector<dcmplx>(2));

    // Calculate transfer matrix elements
    tMatrix[0][0] = 1.0 / t;  // First row, first column
    tMatrix[0][1] = r / t;    // First row, second column
    tMatrix[1][0] = r / t;    // Second row, first column
    tMatrix[1][1] = 1.0 / t;  // Second row, second column

    return tMatrix;
}

/**
 * @details This function calculates the propagation matrix for a specified layer index, 
 * given a transverse wavenumber. The propagation matrix describes the phase 
 * evolution of waves as they propagate through the layer.
 */
std::vector<std::vector<dcmplx>> LayeredMediaUtils::propagationMatrix(int layerIndex, dcmplx krho) const {
    // Validate layer index
    if (layerIndex < 0 || layerIndex >= epsilon.size()) {
        throw std::out_of_range("Invalid layer index for propagation matrix calculation.");
    }

    // Calculate the propagation constant in the layer
    dcmplx kZ = k_L[layerIndex];
    kZ = kZ * kZ - krho * krho;
    kZ = complexSqrt(kZ);
    
    // Prepare the propagation matrix
    std::vector<std::vector<dcmplx>> pMatrix;
    pMatrix.resize(2, std::vector<dcmplx>(2));

    // Calculate propagation matrix elements
    pMatrix[0][0] = exp(-I * kZ * thickness[layerIndex - 1]);
    pMatrix[0][1] = dcmplx(0., 0.);
    pMatrix[1][0] = dcmplx(0., 0.);
    pMatrix[1][1] = exp(I * kZ * thickness[layerIndex - 1]);

    if (std::real(pMatrix[0][0]) > 1e30 || std::imag(pMatrix[0][0]) > 1e30)  pMatrix[0][0] = 1e30;
    if (std::real(pMatrix[0][1]) > 1e30 || std::imag(pMatrix[0][1]) > 1e30)  pMatrix[0][1] = 1e30;
    if (std::real(pMatrix[1][0]) > 1e30 || std::imag(pMatrix[1][0]) > 1e30)  pMatrix[1][0] = 1e30;
    if (std::real(pMatrix[1][1]) > 1e30 || std::imag(pMatrix[1][1]) > 1e30)  pMatrix[1][1] = 1e30;    

    return pMatrix;
}

/**
 * @details This function calculates the cumulative transfer matrix by iteratively multiplying 
 * the propagation matrices and interface transfer matrices for all layers in the structure.
 */
std::vector<std::vector<dcmplx>> LayeredMediaUtils::totalTransferMatrix(dcmplx krho,
                                                                        const std::string& waves) const {

    // Validate wave type
    if (waves != "TE" && waves != "TM") {
        throw std::invalid_argument("Invalid wave type. Use 'TE' or 'TM'.");
    }
    
    std::vector<std::vector<dcmplx>> mTot = this->InterfaceTransferMatrix(0, 1, krho, waves);    
    for (int i = 1; i < zValsInterfaces.size(); i++) {
      std::vector<std::vector<dcmplx>> propMatrix = propagationMatrix(i, krho);
      std::vector<std::vector<dcmplx>> interfaceMatrix = InterfaceTransferMatrix(i, i + 1, krho, waves);
      std::vector<std::vector<dcmplx>> temp = matrixMultiply2D(propMatrix, interfaceMatrix);
      mTot = matrixMultiply2D(mTot, temp);
    }
    return mTot;
}

/**
 * @details This function calculates the overall Fresnel reflection and transmission coefficients 
 * by using the total transfer matrix of the multilayered medium, based on the given 
 * transverse wavenumber and wave polarization.
 */
std::vector<dcmplx> LayeredMediaUtils::totalFresnelCoeffs(dcmplx krho, const std::string& waves, 
                                                                       const std::string& direction) const {
  // Validate wave type
  if (waves != "TE" && waves != "TM") {
    throw std::invalid_argument("Invalid wave type. Use 'TE' or 'TM'.");
  }

  // Validate direction type
  if (direction != "up" && direction != "down") {
    throw std::invalid_argument("Invalid direction type. Use 'up' or 'down'.");
  }

  // Calculate the total transfer matrix
  std::vector<std::vector<dcmplx>> mTot = this->totalTransferMatrix(krho, waves);

  dcmplx totalFresnelReflCoeff(0.0, 0.0);
  dcmplx totalFresnelTransCoeff(0.0, 0.0);
  
  // Extract the reflection and transmission coefficients
  if (direction == "up") {
    // Upgoing wave Prof. Hohenester's book Eq. (8.34)
    totalFresnelReflCoeff  = mTot[1][0] / mTot[0][0];
    totalFresnelTransCoeff = 1. / mTot[0][0];
  } else if (direction == "down") {
    // Downgoing wave Prof. Hohenester's book Eq. (8.36)
    totalFresnelReflCoeff  = -mTot[0][1] / mTot[0][0];
    totalFresnelTransCoeff =  mTot[1][1] + mTot[1][0] * totalFresnelReflCoeff;
  } else {
    throw std::invalid_argument("Invalid direction. Use 'up' or 'down'.");
  }
  
  std::vector<dcmplx> totalFresnelCoeffs = {totalFresnelReflCoeff, totalFresnelTransCoeff};
  
  return totalFresnelCoeffs;
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

/**
 * @details This function calculates the secondary field contributions in a layered medium
 * due to a source located in a specific layer, using the transverse wavenumber
 * and the selected wave polarization (TE or TM).
 */
blitz::Array<dcmplx, 2> LayeredMediaUtils::Secondary(dcmplx krho, int observationLayer, 
                                                     int sourceLayer, const std::string& waves) {

  // Validate wave type
  if (waves != "TE" && waves != "TM") {
    throw std::invalid_argument("Invalid wave type. Use 'TE' or 'TM'.");
  }

  int interfacesNum = this->zValsInterfaces.size();  // Number of interfaces

  // Transfer matrices for each interface
  std::vector<std::vector<std::vector<dcmplx>>> tMatrices(interfacesNum);
  for (int i = 0; i < interfacesNum; ++i) {
    tMatrices[i] = this->InterfaceTransferMatrix(i, i + 1, krho, waves);
  }

  // Compute kZ (propagation constant in each layer)
  std::vector<dcmplx> kZ(interfacesNum + 1);
  for (int i = 0; i < interfacesNum + 1; ++i) {
    dcmplx tmp = k_L[i];
    kZ[i] = complexSqrt(tmp * tmp - krho * krho);
  }

  // Create a 2 x (interfacesNum + 1) matrix
  std::vector<int> range(interfacesNum + 1);
  for (int i = 0; i <= interfacesNum; ++i) {
    range[i] = i;
  }
  std::vector<std::vector<int>> ind(2, range);

  // Interleave the values and remove first and last
  std::vector<int> interleaved;
  for (size_t j = 0; j < ind[0].size(); ++j) {
    interleaved.push_back(ind[0][j]); // Add value from the first row
    interleaved.push_back(ind[1][j]); // Add value from the second row
  }
  interleaved.erase(interleaved.begin());
  interleaved.erase(interleaved.end() - 1);

  int size = 2 * interfacesNum + 2;
  std::vector<bool> ind1(size, true); // Initialize all elements to false
  std::vector<bool> ind2(size, true); // Initialize all elements to false
  ind1[0] = false;                    // Set the first index to true
  ind1[size - 1] = false;             // Set the last index to true
  ind2[1] = false;                    // Set the second index to true
  ind2[size - 2] = false;             // Set the second-to-last index to true

  // Set up matrices for secondary wave coefficients
  blitz::Array<dcmplx, 2> lhs(2 * interfacesNum, 2 * interfacesNum + 2, blitz::ColumnMajorArray<2>());
  blitz::Array<dcmplx, 2> rhs(2 * interfacesNum, 2 * interfacesNum + 2, blitz::ColumnMajorArray<2>());

  // Initialize propagation constants
  std::vector<dcmplx> fac(interfacesNum);
  fac[0] = 1.0;  // First element is 0
  for (int i = 1; i < interfacesNum; ++i) {
    fac[i] = exp(I * kZ[i] * thickness[i - 1]);
  }
  fac.push_back(1.0);  // Add the last element as 0

  for (int i = 0; i < interfacesNum; ++i) {
    // Extract transfer matrix
    auto M = tMatrices[i];

    // Define sub-matrix indices
    std::vector<int> k1 = {2 * i, 2 * i + 1};
    // Compute k2 indices (zero-based)
    std::vector<int> k2 = {2 * i + 2, 2 * i + 3};

    // Fill lhs and rhs matrices using calculated indices
    lhs(k1[0], k1[0]) = fac[i];
    lhs(k1[0], k1[1]) = 0.0;
    lhs(k1[1], k1[0]) = 0.0;
    lhs(k1[1], k1[1]) = 1.0;

    rhs(k1[0], k1[0]) = 1.0;
    rhs(k1[0], k1[1]) = 0.0;
    rhs(k1[1], k1[0]) = 0.0;
    rhs(k1[1], k1[1]) = 0.0;

    lhs(k1[0], k2[0]) = -M[0][0];
    lhs(k1[0], k2[1]) = -M[0][1] * fac[i + 1];
    lhs(k1[1], k2[0]) = -M[1][0];
    lhs(k1[1], k2[1]) = -M[1][1] * fac[i + 1];

    rhs(k1[0], k2[0]) = -0.0;
    rhs(k1[0], k2[1]) = -M[0][1];
    rhs(k1[1], k2[0]) = -0.0;
    rhs(k1[1], k2[1]) = -M[1][1];
  }

  int lhsCols = std::count(ind1.begin(), ind1.end(), true); // Count how many 'true' values in ind1
  int rhsCols = std::count(ind2.begin(), ind2.end(), true); // Count how many 'true' values in ind2
  blitz::Array<dcmplx, 2> lhsNew(lhs.extent(0), lhsCols, blitz::ColumnMajorArray<2>());
  blitz::Array<dcmplx, 2> rhsNew(rhs.extent(0), rhsCols, blitz::ColumnMajorArray<2>());
  blitz::Array<dcmplx, 2> sysSol(lhs.extent(0), rhsCols, blitz::ColumnMajorArray<2>());

  int colIndex1 = 0;
  for (int i = 0; i < ind1.size(); ++i) {
    if (ind1[i]) { // If the value is true
        lhsNew(blitz::Range::all(), colIndex1++) = lhs(blitz::Range::all(), i);
    }
  }

  int colIndex2 = 0;
  for (int i = 0; i < ind2.size(); ++i) {
    if (ind2[i]) { // If the value is true
        rhsNew(blitz::Range::all(), colIndex2++) = rhs(blitz::Range::all(), i);
    }
  }

  // Solve for the secondary fields
  int arrayN = lhsNew.extent(0);
  blitz::Array<int, 1> iPiv(arrayN);
  int info = 0;
  char trans = 'N';
  blitz::Array<dcmplx, 2> rightSides(arrayN, rhsNew.extent(1), blitz::ColumnMajorArray<2>());
  
  rightSides = rhsNew; // Copy the rhs to rightSides for solving

  // Perform LU factorization
  zgetrf_(&arrayN, &arrayN, &(lhsNew(0, 0)), &arrayN, &(iPiv(0)), &info);
  
  // Solve the system
  zgetrs_(&trans, &arrayN, &arrayN, &(lhsNew(0, 0)), &arrayN, &(iPiv(0)), &(rightSides(0, 0)), &arrayN, &info);

  // Extract solution from rightSides back to the result
  for (int i = 0; i < rhsNew.extent(1); i++) {
    sysSol(blitz::Range::all(), i) = -rightSides(blitz::Range::all(), i);
  }

  // Select the rows based on interleaved
  std::vector<int> rowIndices;
  for (int i = 0; i < interleaved.size(); ++i) {
    if (interleaved[i] == observationLayer) {
        rowIndices.push_back(i); // Add the row index if the value is true
    }
  }

  // Select the columns based on interleaved
  std::vector<int> colIndices;
  for (int i = 0; i < interleaved.size(); ++i) {
    if (interleaved[i] == sourceLayer) {
        colIndices.push_back(i); // Add the column index if the value is true
    }
  }

  // Define the result sub-matrix
  blitz::Array<dcmplx, 2> result(rowIndices.size(), colIndices.size(), blitz::ColumnMajorArray<2>());

  // Extract the submatrix from x using row and column indices
  int rowIndex = 0;
  for (int i : rowIndices) {
    int colIndex = 0;
    for (int j : colIndices) {
        result(rowIndex, colIndex++) = sysSol(i, j);
    }
    ++rowIndex;
  }

#ifdef DEBUG
  for (int i = 0; i < result.extent(0); ++i) { // Rows
    for (int j = 0; j < result.extent(1); ++j) { // Columns
        std::cout << "(" << std::real(result(i, j)) << " + "
                  << std::imag(result(i, j)) << "i) ";
    }
    std::cout << std::endl; // Newline after each row
  }
#endif
  
  return result;
}

/**
 * @details This function calculates the field values at specified observation points
 * in a layered medium, given a source position, transverse wavenumber, and
 * wave polarization (TE or TM).
 */
std::vector<dcmplx> LayeredMediaUtils::Fields(dcmplx krho, std::vector<double> observationPointsZ, 
                                              double sourcePointZ, const std::string& waves) {
  // Resulting electric field (TE or TM)
  std::vector<dcmplx> fields(observationPointsZ.size(), dcmplx(0.0, 0.0));

  int interfacesNum = this->zValsInterfaces.size();  // Number of interfaces

  // Get layer indices for observation points and source point
  std::vector<int> obsPointsLayersIndices(observationPointsZ.size(), -1); // Indices for observation points
  int srcPointLayerIndex = GetLayerIndex(rvec(0., 0., sourcePointZ));     // Index for the source point
  for (int i = 0; i < observationPointsZ.size(); i++){
    obsPointsLayersIndices[i] = GetLayerIndex(rvec(0., 0., observationPointsZ[i]));
  }
  
  // Wave propagation direction
  const int dir[] = {1, -1};

  // Loop over unique indices of observation layers
  std::unordered_set<int> temp(obsPointsLayersIndices.begin(), obsPointsLayersIndices.end());
  std::vector<int> uniqueLayersIndices(temp.begin(), temp.end());
  std::sort(uniqueLayersIndices.begin(), uniqueLayersIndices.end());
  temp.clear();
  
  for (int it : uniqueLayersIndices) {
    // Coefficients for secondary waves
    auto waveSec = Secondary(krho, it, srcPointLayerIndex, waves);
    
    // Wavenumbers' z components
    dcmplx kzLayerObs = k_L[it];
    dcmplx kzLayerSrc = k_L[srcPointLayerIndex];
    kzLayerObs = kzLayerObs * kzLayerObs - krho * krho; kzLayerObs = complexSqrt(kzLayerObs);
    kzLayerSrc = kzLayerSrc * kzLayerSrc - krho * krho; kzLayerSrc = complexSqrt(kzLayerSrc);

    // Indices corresponding to the current layer
    std::vector<bool> isInCurrentLayer(observationPointsZ.size(), false);

    for (int j = 0; j < obsPointsLayersIndices.size(); ++j) {
      if (obsPointsLayersIndices[j] == it) {
          isInCurrentLayer[j] = true;
      }
    }

    // Loop over wave directions
    for (int dObs = 0; dObs < waveSec.extent(0); dObs++) {
      for (int dSrc = 0; dSrc < waveSec.extent(1); dSrc++) {
        int dirObs = dir[dObs]; if (it == 0) dirObs = -1;
        int dirSrc = dir[dSrc]; if (srcPointLayerIndex == interfacesNum + 1) dirSrc = -1;

        // Compute distances to interfaces
        std::vector<double> distObs;
        double distSrc = 0.0;
        for (int i = 0; i < observationPointsZ.size(); i++) {
          if (isInCurrentLayer[i]){
            if (dirObs == 1) {
              // Distance from the lower interface
              distObs.push_back(observationPointsZ[i] - zValsInterfaces[it - 1]);
            }
            else {
              // Distance from the upper interface
              distObs.push_back(zValsInterfaces[it] - observationPointsZ[i]);
            }
          }
        }
        
        if (dirSrc == -1) {
          // Distance from the lower interface
          distSrc = sourcePointZ - zValsInterfaces[srcPointLayerIndex - 1];
        }
        else {
          // Distance from the upper interface
          distSrc = zValsInterfaces[srcPointLayerIndex] - sourcePointZ;
        }

        // Accumulate secondary fields
        int counter = 0;
        for (int i = 0; i < observationPointsZ.size(); ++i) {
          if (isInCurrentLayer[i]) {
            fields[i] += waveSec(dObs, dSrc) * exp(I * (kzLayerObs * distObs[counter] + kzLayerSrc * distSrc));
            counter++;
          }
        }
      }
    }

    // Add primary field if needed
    if (it == srcPointLayerIndex) {
      for (int i = 0; i < observationPointsZ.size(); ++i) {
        if (isInCurrentLayer[i]) {
          fields[i] += exp(I * kzLayerObs * std::abs(observationPointsZ[i] - sourcePointZ));
        }
      }
    }
  }

#ifdef DEBUG
  for (const auto& elem : fields){
    std::cout << elem << ' ';
    std::cout << '\n';
  }
#endif

  return fields;
}

/**
 * @details This function multiplies two 2D complex matrices and returns the resulting matrix.
 * It ensures that the input matrices have compatible dimensions for multiplication.
 */
std::vector<std::vector<dcmplx>> LayeredMediaUtils::matrixMultiply2D(const std::vector<std::vector<dcmplx>>& x,
                                                                     const std::vector<std::vector<dcmplx>>& y) const{

  if (x[0].size() != y.size()) {
      throw std::invalid_argument("Matrix dimensions do not match for multiplication.");
  }

  size_t rows = x.size();
  size_t cols = y[0].size();
  size_t common = y.size();

  std::vector<std::vector<dcmplx>> z(rows, std::vector<dcmplx>(cols, 0.0));

  for (size_t i = 0; i < rows; ++i) {
      for (size_t j = 0; j < cols; ++j) {
          for (size_t k = 0; k < common; ++k) {
              z[i][j] += x[i][k] * y[k][j];
          }
      }
  }
  return z;
}

/**
 * @details This function multiplies corresponding elements of two 3D complex matrices,
 * assuming a 2x2 structure across multiple depth slices.
 */
std::vector<std::vector<std::vector<dcmplx>>> LayeredMediaUtils::matrixMultiply3D(
                                                        const std::vector<std::vector<std::vector<dcmplx>>>& x,
                                                        const std::vector<std::vector<std::vector<dcmplx>>>& y) const{

  size_t depth = x[0][0].size(); // Number of slices
  std::vector<std::vector<std::vector<dcmplx>>> z(2, std::vector<std::vector<dcmplx>>(2, std::vector<dcmplx>(depth, 0.0)));

  for (size_t d = 0; d < depth; ++d) {
      z[0][0][d] = x[0][0][d] * y[0][0][d] + x[0][1][d] * y[1][0][d];
      z[1][0][d] = x[1][0][d] * y[0][0][d] + x[1][1][d] * y[1][0][d];
      z[0][1][d] = x[0][0][d] * y[0][1][d] + x[0][1][d] * y[1][1][d];
      z[1][1][d] = x[1][0][d] * y[0][1][d] + x[1][1][d] * y[1][1][d];
  }

  return z;
}

/**
 * @details This function determines whether the input matrices are 2D or 3D and calls 
 * the appropriate multiplication function (`matrixMultiply2D` or `matrixMultiply3D`). 
 * The results are printed to the console.
 */
void LayeredMediaUtils::matrixMultiplication(const std::vector<std::vector<std::vector<dcmplx>>>& x,
                                             const std::vector<std::vector<std::vector<dcmplx>>>& y) const{

  if (x[0][0].size() == 1) { // 2D Case
      std::vector<std::vector<dcmplx>> x2D = {{x[0][0][0], x[0][1][0]}, {x[1][0][0], x[1][1][0]}};
      std::vector<std::vector<dcmplx>> y2D = {{y[0][0][0], y[0][1][0]}, {y[1][0][0], y[1][1][0]}};

      std::vector<std::vector<dcmplx>> result = matrixMultiply2D(x2D, y2D);

      std::cout << "2D Multiplication Result:\n";
      for (const auto& row : result) {
          for (dcmplx val : row) {
              std::cout << val << " ";
          }
          std::cout << std::endl;
      }
  } else { // Higher-dimensional case
      std::vector<std::vector<std::vector<dcmplx>>> result = matrixMultiply3D(x, y);

      std::cout << "3D Multiplication Result:\n";
      for (size_t d = 0; d < result[0][0].size(); ++d) {
          std::cout << "Slice " << d + 1 << ":\n";
          for (size_t i = 0; i < 2; ++i) {
              for (size_t j = 0; j < 2; ++j) {
                  std::cout << result[i][j][d] << " ";
              }
              std::cout << std::endl;
          }
          std::cout << std::endl;
      }
  }
}

/**
 * @details Converts Cartesian coordinates to layered medium coordinates for different 
 * interaction types regarding the calculation of the quasistatic part of Green's function.
 */
layeredCoords LayeredMediaUtils::cartesianToLayeredQS(const rvec &r1, const rvec &r2, 
                                                      const std::string& type) const
{
    layeredCoords out{};
    out.r  = 0.0; // Radial distance in xy‐plane
    out.Z1 = 0.0; // Vertical distance for lower interface (intraIn) or Z for intraTB/inter
    out.Z2 = 0.0; // Vertical distance for upper interface (intraIn)
    out.R1 = 0.0; // Total distance for lower interface (intraIn) or R for intraTB/inter
    out.R2 = 0.0; // Total distance for upper interface

    // Radial distance in xy‐plane:
    double dx = r1[0] - r2[0];
    double dy = r1[1] - r2[1];
    out.r = std::sqrt(dx * dx + dy * dy);
    
    // Vertical and total distances:
    if (type == "intraTB") {
        // For intra-layer interaction when layer is the top/bottom layer
        int indexLayer = GetLayerIndex(r1);
        if (indexLayer < 0 || indexLayer > epsilon.size() - 1) {
            throw std::out_of_range("Invalid layer index for intraTB interaction.");
        }
        double interfaceZ = (indexLayer == 0 ? 
                             zValsInterfaces.front() : zValsInterfaces.back());
        out.Z1 = std::fabs(r1[2] + r2[2] - 2.0 * interfaceZ);
        out.R1 = std::sqrt(out.r * out.r + out.Z1 * out.Z1);
    }
    else if (type == "intraIn") {
        // For intra-layer interaction when layer is NOT the top/bottom layer
        int indexLayer = GetLayerIndex(r1);
        if (indexLayer <= 0 || indexLayer >= epsilon.size() - 1) {
            throw std::out_of_range("Invalid layer index for intraIn interaction.");
        }
        double lowerInterfaceZ = zValsInterfaces[indexLayer - 1];
        double upperInterfaceZ = zValsInterfaces[indexLayer];
        out.Z1 = std::fabs(r1[2] + r2[2] - 2.0 * lowerInterfaceZ);
        out.Z2 = std::fabs(r1[2] + r2[2] - 2.0 * upperInterfaceZ);
        out.R1 = std::sqrt(out.r * out.r + out.Z1 * out.Z1);
        out.R2 = std::sqrt(out.r * out.r + out.Z2 * out.Z2);
    }
    else if (type == "inter") {
        // For inter-layer interaction (source & observation layers are different)
        out.Z1 = std::fabs(r1[2] - r2[2]);
        out.R1 = std::sqrt(out.r * out.r + out.Z1 * out.Z1);
    }
    else {
        throw std::invalid_argument("Invalid interaction type. Use 'intraTB', 'intraIn', or 'inter'.");
    }
    if (out.r < 1.0e-10) out.r = 1.0e-10; // Avoid division by zero in further calculations
    return out;
}

/**
 * @details Converts Cartesian coordinates to layered medium coordinates for different 
 * interaction types regarding the calculation of the smooth part of Green's function.
 */
layeredCoords LayeredMediaUtils::cartesianToLayeredTab(const rvec &r1, const rvec &r2, 
                                                       const std::string& type) const
{
    layeredCoords out{};
    out.r  = 0.0; // Radial distance in xy‐plane
    out.Z1 = 0.0; // Vertical distance for lower interface (intraIn) or Z for intraTB/inter
    out.Z2 = 0.0; // Vertical distance for upper interface (intraIn)
    out.R1 = 0.0; // Total distance for lower interface (intraIn) or R for intraTB/inter
    out.R2 = 0.0; // Total distance for upper interface

    // Radial distance in xy‐plane:
    double dx = r1[0] - r2[0];
    double dy = r1[1] - r2[1];
    out.r = std::sqrt(dx * dx + dy * dy);
    
    // Vertical and total distances:
    if (type == "intraTB") { // For intra-layer interaction when layer is the top/bottom layer
        
        // Check if the layer is the top or bottom layer
        int indexLayer = GetLayerIndex(r1);
        if (indexLayer < 0 || indexLayer > epsilon.size() - 1) {
            throw std::out_of_range("Invalid layer index for intraTB interaction.");
        }

        // Calculate the vertical distance based on the layer index
        double interfaceZ = (indexLayer == 0 ? 
                             zValsInterfaces.front() : zValsInterfaces.back());
        out.Z1 = std::fabs(r1[2] + r2[2] - 2.0 * interfaceZ);
    }
    else if (type == "intraIn") { // For intra-layer interaction (layer is NOT the first/last)

        // Check if the layer is a layer of the interior
        int indexLayer = GetLayerIndex(r1);
        if (indexLayer <= 0 || indexLayer >= epsilon.size() - 1) {
            throw std::out_of_range("Invalid layer index for intraIn interaction.");
        }
        
        // Calculate the vertical distance based on the layer index
        double lowerInterfaceZ = zValsInterfaces[indexLayer - 1];
        double upperInterfaceZ = zValsInterfaces[indexLayer];
        out.Z1 = r1[2] - r2[2] + upperInterfaceZ - lowerInterfaceZ;
        out.Z2 = r1[2] + r2[2] - 2.0 * lowerInterfaceZ;
    }
    else if (type == "inter") { // For inter-layer interaction (source & observation layers are different)
        out.Z1 = r1[2]; out.Z2 = r2[2];
    }
    else {
        throw std::invalid_argument("Invalid interaction type. Use 'intraTB', 'intraIn', or 'inter'.");
    }
    return out;
}