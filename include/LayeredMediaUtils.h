/**
 * @file LayeredMediaUtils.h
 * @brief This header declares utility structures and algorithms used throughout the
 * layered-media scattering solver. The utilities cover:
 * - Indexing of layers from Cartesian z (GetLayerIndex)
 * - Fresnel coefficients, transfer and propagation matrices for TE/TM (FresnelCoeff,
 * InterfaceTransferMatrix, propagationMatrix, totalTransferMatrix)
 * - Total reflection/transmission of a multilayer stack (totalFresnelCoeffs)
 * - Source/observation coordinate transforms for the fast-tabulated Green's
 * functions (cartesianToLayeredTab, cartesianToLayeredQS)
 * - Discrete range slicing and gridding for tabulation (slice, computeGrid)
 * - Reading helper utilities for test data (readPositionsFromFile)
 * See the .cpp for algorithmic details and numerical conventions.
 */

#ifndef LAYERED_MEDIA_UTILS_H
#define LAYERED_MEDIA_UTILS_H

#include <map>
#include <stdexcept>         // For exception handling
#include <functional>        // For std::function
#include <vector>
#include "globals.h"
#include <fstream>

// ---------------------------------------------------------------------------
//                     Layered media utility structures
// ---------------------------------------------------------------------------
/**
 * @struct PositionGroup
 * @brief Structure representing a group of spatial positions categorized by layer indices.
 */
struct PositionGroup {
  size_t i1;                ///< Layer index1
  size_t i2;                ///< Layer index2
  std::vector<rvec> pos1;   ///< Grouped positions
  std::vector<rvec> pos2;   ///< Grouped positions
  std::vector<size_t> ind1; ///< Indices of grouped positions
  std::vector<size_t> ind2; ///< Indices of grouped positions
};

/**
 * @struct Range
 * @brief Structure representing the spatial ranges in a layered medium tabulation.
 */
struct Range {
  int i1;                 ///< Layer index1
  int i2;                 ///< Layer index2
  std::vector<double> r;  ///< Radial range
  std::vector<double> z1; ///< Z1 range
  std::vector<double> z2; ///< Z2 range (optional)
};

/**
 * @struct Grid
 * @brief Structure representing the grid points used for tabulation in a layered medium.
 */
struct Grid {
  int i1;                 ///< Observation layer index
  int i2;                 ///< Source layer index
  std::vector<double> r;  ///< Radial grid points
  std::vector<double> z1; ///< Z1 grid points
  std::vector<double> z2; ///< Z2 grid points (optional)
};

/**
 * @struct layeredCoords
 * @brief Structure representing stratified coordinates in a layered medium.
 */
struct layeredCoords {
    double r{0.0};  ///< Radial distance in xy-plane
    double Z1{0.0}; ///< Vertical distance for lower interface (intraIn) or Z for intraTB/inter
    double Z2{0.0}; ///< Vertical distance for upper interface (intraIn)
    double R1{0.0}; ///< Total distance for lower interface (intraIn) or R for intraTB/inter
    double R2{0.0}; ///< Total distance for upper interface
};

/**
 * @struct Data2D
 * @brief Structure representing the 2D spatial data used for tabulation in a layered medium.
 */
struct Data2D {
  std::vector<double> rVec;                     ///< Radial coordinate sampling points.
  std::vector<std::vector<double>> rMatrix;     ///< 2D grid of radial coordinates.
  std::vector<double> z1Vec;                    ///< Z1 coordinate sampling points.
  std::vector<std::vector<double>> z1Matrix;    ///< 2D grid of Z1 coordinates.
  std::vector<double> z2Vec;                    ///< Z2 coordinate sampling points.
  std::vector<std::vector<double>> z2Matrix;    ///< 2D grid of Z2 coordinates.
  std::vector<int> indVec;                      ///< Flat index mapping for grid access.
  std::vector<std::complex<int>> indVecComplex; ///< Complex-valued index mapping.
};

/**
 * @struct Data3D
 * @brief Structure representing the 3D spatial data used for tabulation in a layered medium.
 */
struct Data3D {
  std::vector<double> rVec;                               ///< Radial coordinate sampling points.
  std::vector<std::vector<std::vector<double>>> rMatrix;  ///< 3D grid of radial coordinates.
  std::vector<double> z1Vec;                              ///< Z1 coordinate sampling points.
  std::vector<std::vector<std::vector<double>>> z1Matrix; ///< 3D grid of Z1 coordinates.
  std::vector<double> z2Vec;                              ///< Z2 coordinate sampling points.
  std::vector<std::vector<std::vector<double>>> z2Matrix; ///< 3D grid of Z2 coordinates.
  std::vector<int> indVec;                                ///< Flat index mapping for grid access.
  std::vector<std::tuple<int, int, int>> indVecTriplet;   ///< Triplet index mapping (rho, Z1, Z2).
};

/**
 * @struct CoeffsIntra
 * @brief Interpolation coefficients of the intra-layer interactions of Green's function.
 */
struct CoeffsIntra {
  std::tuple<std::vector<dcmplx>, std::vector<int>, 
             std::vector<int>, std::vector<double>, 
             std::vector<double>> te;               ///< TE intra-layer coefficients.
  std::tuple<std::vector<dcmplx>, std::vector<int>,
             std::vector<int>, std::vector<double>, 
             std::vector<double>> tm;               ///< TM intra-layer coefficients.
  std::tuple<std::vector<dcmplx>, std::vector<int>, 
             std::vector<int>, std::vector<double>, 
             std::vector<double>> tez;              ///< TE–Z intra-layer coefficients.
  std::tuple<std::vector<dcmplx>, std::vector<int>, 
             std::vector<int>, std::vector<double>, 
             std::vector<double>> tmz;              ///< TM–Z intra-layer coefficients.
  std::tuple<std::vector<dcmplx>, std::vector<int>, 
             std::vector<int>, std::vector<double>, 
             std::vector<double>> tezz;             ///< TE–ZZ intra-layer coefficients.
  std::tuple<std::vector<dcmplx>, std::vector<int>, 
             std::vector<int>, std::vector<double>, 
             std::vector<double>> tmzz;             ///< TM–ZZ intra-layer coefficients.
  std::tuple<std::vector<dcmplx>, std::vector<int>, 
             std::vector<int>, std::vector<double>, 
             std::vector<double>> tes;              ///< TE–S intra-layer coefficients.
  std::tuple<std::vector<dcmplx>, std::vector<int>, 
             std::vector<int>, std::vector<double>, 
             std::vector<double>> tms;              ///< TM–S intra-layer coefficients.
  std::tuple<std::vector<dcmplx>, std::vector<int>, 
             std::vector<int>, std::vector<double>, 
             std::vector<double>> ter;              ///< TE–R intra-layer coefficients.
  std::tuple<std::vector<dcmplx>, std::vector<int>, 
             std::vector<int>, std::vector<double>, 
             std::vector<double>> tmr;              ///< TM–R intra-layer coefficients.
  std::tuple<std::vector<dcmplx>, std::vector<int>, 
             std::vector<int>, std::vector<double>, 
             std::vector<double>> terz;             ///< TE–RZ intra-layer coefficients.
  std::tuple<std::vector<dcmplx>, std::vector<int>, 
             std::vector<int>, std::vector<double>, 
             std::vector<double>> tmrz;             ///< TM–RZ intra-layer coefficients.
};

/**
 * @struct CoeffsInter
 * @brief Interpolation coefficients of the inter-layer interactions of Green's function.
 */
struct CoeffsInter {
  std::tuple<std::vector<double>, std::vector<double>, 
             std::vector<double>, std::vector<dcmplx>> te;    ///< TE inter-layer coefficients.
  std::tuple<std::vector<double>, std::vector<double>, 
             std::vector<double>, std::vector<dcmplx>> tm;    ///< TM inter-layer coefficients.
  std::tuple<std::vector<double>, std::vector<double>,  
             std::vector<double>, std::vector<dcmplx>> tez1;  ///< TE–Z1 inter-layer coefficients.
  std::tuple<std::vector<double>, std::vector<double>,  
             std::vector<double>, std::vector<dcmplx>> tmz1;  ///< TM–Z1 inter-layer coefficients.
  std::tuple<std::vector<double>, std::vector<double>,  
             std::vector<double>, std::vector<dcmplx>> tez2;  ///< TE–Z2 inter-layer coefficients.
  std::tuple<std::vector<double>, std::vector<double>,  
             std::vector<double>, std::vector<dcmplx>> tmz2;  ///< TM–Z2 inter-layer coefficients.
  std::tuple<std::vector<double>, std::vector<double>,  
             std::vector<double>, std::vector<dcmplx>> tezz;  ///< TE–ZZ inter-layer coefficients.
  std::tuple<std::vector<double>, std::vector<double>,  
             std::vector<double>, std::vector<dcmplx>> tmzz;  ///< TM–ZZ inter-layer coefficients.
  std::tuple<std::vector<double>, std::vector<double>,  
             std::vector<double>, std::vector<dcmplx>> tes;   ///< TE–S inter-layer coefficients.
  std::tuple<std::vector<double>, std::vector<double>,  
             std::vector<double>, std::vector<dcmplx>> tms;   ///< TM–S inter-layer coefficients.
  std::tuple<std::vector<double>, std::vector<double>,  
             std::vector<double>, std::vector<dcmplx>> ter;   ///< TE–R inter-layer coefficients.
  std::tuple<std::vector<double>, std::vector<double>,  
             std::vector<double>, std::vector<dcmplx>> tmr;   ///< TM–R inter-layer coefficients.
  std::tuple<std::vector<double>, std::vector<double>,  
             std::vector<double>, std::vector<dcmplx>> terz1; ///< TE–RZ1 inter-layer coefficients.
  std::tuple<std::vector<double>, std::vector<double>,  
             std::vector<double>, std::vector<dcmplx>> tmrz1; ///< TM–RZ1 inter-layer coefficients.
  std::tuple<std::vector<double>, std::vector<double>,  
             std::vector<double>, std::vector<dcmplx>> terz2; ///< TE–RZ2 inter-layer coefficients.
  std::tuple<std::vector<double>, std::vector<double>,  
             std::vector<double>, std::vector<dcmplx>> tmrz2; ///< TM–RZ2 inter-layer coefficients.
  
  std::vector<double> scales;                                 ///< Scaling factors.
};

/**
 * @struct IntegrandIntra
 * @brief Structure representing the integrand of the intra-layer interactions of Green's function.
 */
struct IntegrandIntra {
  std::vector<std::vector<dcmplx>> te;    ///< TE intra-layer integrand.
  std::vector<std::vector<dcmplx>> tm;    ///< TM intra-layer integrand.
  std::vector<std::vector<dcmplx>> tez;   ///< TE–Z intra-layer integrand.
  std::vector<std::vector<dcmplx>> tmz;   ///< TM–Z intra-layer integrand.
  std::vector<std::vector<dcmplx>> tezz;  ///< TE–ZZ intra-layer integrand.
  std::vector<std::vector<dcmplx>> tmzz;  ///< TM–ZZ intra-layer integrand.
  std::vector<std::vector<dcmplx>> tes;   ///< TE–S intra-layer integrand.
  std::vector<std::vector<dcmplx>> tms;   ///< TM–S intra-layer integrand.
  std::vector<std::vector<dcmplx>> ter;   ///< TE–R intra-layer integrand.
  std::vector<std::vector<dcmplx>> tmr;   ///< TM–R intra-layer integrand.
  std::vector<std::vector<dcmplx>> terz;  ///< TE–RZ intra-layer integrand.
  std::vector<std::vector<dcmplx>> tmrz;  ///< TM–RZ intra-layer integrand.

  /**
  * @brief Overloaded subtraction operator.
  * @param other The other IntegrandIntra to subtract.
  * @return The resulting IntegrandIntra after subtraction.
  */
  IntegrandIntra operator-(const IntegrandIntra& other) const {
    IntegrandIntra result = *this;
    subtractVector(result.te,   other.te  );
    subtractVector(result.tm,   other.tm  );
    subtractVector(result.tez,  other.tez );
    subtractVector(result.tmz,  other.tmz );
    subtractVector(result.tezz, other.tezz);
    subtractVector(result.tmzz, other.tmzz);
    subtractVector(result.tes,  other.tes );
    subtractVector(result.tms,  other.tms );
    subtractVector(result.ter,  other.ter );
    subtractVector(result.tmr,  other.tmr );
    subtractVector(result.terz, other.terz);
    subtractVector(result.tmrz, other.tmrz);
    return result;
  }

  /**
   * @brief Overloaded multiplication operator (scalar * structure).
   * @param scalar The scalar to multiply.
   * @param obj The IntegrandIntra to multiply.
   * @return The resulting IntegrandIntra after multiplication.
   */
  friend IntegrandIntra operator*(const dcmplx& scalar, const IntegrandIntra& obj) {
      IntegrandIntra result = obj;
      multiplyVector(result.te,   scalar);
      multiplyVector(result.tm,   scalar);
      multiplyVector(result.tez,  scalar);
      multiplyVector(result.tmz,  scalar);
      multiplyVector(result.tezz, scalar);
      multiplyVector(result.tmzz, scalar);
      multiplyVector(result.tes,  scalar);
      multiplyVector(result.tms,  scalar);
      multiplyVector(result.ter,  scalar);
      multiplyVector(result.tmr,  scalar);
      multiplyVector(result.terz, scalar);
      multiplyVector(result.tmrz, scalar);
      return result;
  }

  /**
   * @brief Function to create a flattened structure.
   * @param indVecComplex The indices of the flattened structure.
   * @return The flattened structure.
   */
  IntegrandIntra flatten(const std::vector<std::complex<int>>& indVecComplex) const {
    IntegrandIntra result;

    // Process each member
    result.te   = flattenVector(te,   indVecComplex);
    result.tm   = flattenVector(tm,   indVecComplex);
    result.tez  = flattenVector(tez,  indVecComplex);
    result.tmz  = flattenVector(tmz,  indVecComplex);
    result.tezz = flattenVector(tezz, indVecComplex);
    result.tmzz = flattenVector(tmzz, indVecComplex);
    result.tes  = flattenVector(tes,  indVecComplex);
    result.tms  = flattenVector(tms,  indVecComplex);
    result.ter  = flattenVector(ter,  indVecComplex);
    result.tmr  = flattenVector(tmr,  indVecComplex);
    result.terz = flattenVector(terz, indVecComplex);
    result.tmrz = flattenVector(tmrz, indVecComplex);

    return result;
}

  private:
  /**
   * @brief Helper function to subtract two 2D vectors element-wise.
   * @param a The first vector.
   * @param b The second vector.
   */
  static void subtractVector(std::vector<std::vector<dcmplx>>& a, 
                             const std::vector<std::vector<dcmplx>>& b) {
    if (a.size() != b.size()) return; // Ensure sizes match
      for (size_t i = 0; i < a.size(); ++i) {
          if (a[i].size() != b[i].size()) continue;
          for (size_t j = 0; j < a[i].size(); ++j) {
              a[i][j] -= b[i][j];
          }
      }
    }

  /**
   * @brief Helper function to multiply a 2D vector by a scalar.
   * @param vec The vector to multiply.
   * @param scalar The scalar to multiply by.
   */
  static void multiplyVector(std::vector<std::vector<dcmplx>>& vec, const dcmplx& scalar) {
    for (size_t i = 0; i < vec.size(); ++i) {
        for (size_t j = 0; j < vec[i].size(); ++j) {
            vec[i][j] *= scalar;
        }
    }
  }

  /**
   * @brief Helper function to extract selected elements from a 2D vector.
   * @param matrix The 2D vector to extract from.
   * @param indVecComplex The indices of the elements to extract.
   */
  std::vector<std::vector<dcmplx>> flattenVector(
      const std::vector<std::vector<dcmplx>>& matrix,
      const std::vector<std::complex<int>>& indVecComplex) const {
      std::vector<std::vector<dcmplx>> flattened;
      for (const auto& coord : indVecComplex) {
        int i = coord.real(); // Extract row index
        int j = coord.imag(); // Extract column index
        flattened.push_back({ getValue(matrix, i, j) }); // Store as 1-column 2D vector
      }
      return flattened;
    }

  /**
   * @brief Helper function to safely access matrix elements.
   * @param matrix The 2D vector to access.
   * @param i The row index.
   * @param j The column index.
   */
  dcmplx getValue(const std::vector<std::vector<dcmplx>>& matrix, int i, int j) const {
    if (i >= 0 && i < static_cast<int>(matrix.size()) && 
        j >= 0 && j < static_cast<int>(matrix[i].size())) {
      return matrix[i][j];
    }
      return dcmplx(0, 0); // Default to zero if out of bounds
  }
};

/**
 * @struct IntegrandInter
 * @brief Structure representing the integrand of the inter-layer interactions of Green's function.
 */
struct IntegrandInter {
  std::vector<std::vector<std::vector<dcmplx>>> te;     ///< TE inter-layer integrand.
  std::vector<std::vector<std::vector<dcmplx>>> tm;     ///< TM inter-layer integrand.
  std::vector<std::vector<std::vector<dcmplx>>> tez1;   ///< TE–Z1 inter-layer integrand.
  std::vector<std::vector<std::vector<dcmplx>>> tmz1;   ///< TM–Z1 inter-layer integrand.
  std::vector<std::vector<std::vector<dcmplx>>> tez2;   ///< TE–Z2 inter-layer integrand.
  std::vector<std::vector<std::vector<dcmplx>>> tmz2;   ///< TM–Z2 inter-layer integrand.
  std::vector<std::vector<std::vector<dcmplx>>> tezz;   ///< TE–ZZ inter-layer integrand.
  std::vector<std::vector<std::vector<dcmplx>>> tmzz;   ///< TM–ZZ inter-layer integrand.
  std::vector<std::vector<std::vector<dcmplx>>> tes;    ///< TE–S inter-layer integrand.
  std::vector<std::vector<std::vector<dcmplx>>> tms;    ///< TM–S inter-layer integrand.
  std::vector<std::vector<std::vector<dcmplx>>> ter;    ///< TE–R inter-layer integrand.
  std::vector<std::vector<std::vector<dcmplx>>> tmr;    ///< TM–R inter-layer integrand.
  std::vector<std::vector<std::vector<dcmplx>>> terz1;  ///< TE–RZ1 inter-layer integrand.
  std::vector<std::vector<std::vector<dcmplx>>> tmrz1;  ///< TM–RZ1 inter-layer integrand.
  std::vector<std::vector<std::vector<dcmplx>>> terz2;  ///< TE–RZ2 inter-layer integrand.
  std::vector<std::vector<std::vector<dcmplx>>> tmrz2;  ///< TM–RZ2 inter-layer integrand.

  /** 
   * @brief Overloaded subtraction operator.
   * @param other The other IntegrandInter to subtract.
   * @return The resulting IntegrandInter after subtraction.
   */
  IntegrandInter operator-(const IntegrandInter& other) const {
    IntegrandInter result = *this;
    subtractVector(result.te,    other.te   );
    subtractVector(result.tm,    other.tm   );
    subtractVector(result.tez1,  other.tez1 );
    subtractVector(result.tmz1,  other.tmz1 );
    subtractVector(result.tez2,  other.tez2 );
    subtractVector(result.tmz2,  other.tmz2 );
    subtractVector(result.tezz,  other.tezz );
    subtractVector(result.tmzz,  other.tmzz );
    subtractVector(result.tes,   other.tes  );
    subtractVector(result.tms,   other.tms  );
    subtractVector(result.ter,   other.ter  );
    subtractVector(result.tmr,   other.tmr  );
    subtractVector(result.terz1, other.terz1);
    subtractVector(result.tmrz1, other.tmrz1);
    subtractVector(result.terz2, other.terz2);
    subtractVector(result.tmrz2, other.tmrz2);
    return result;
  }

  /**
   * @brief Overloaded multiplication operator (scalar * structure).
   * @param scalar The scalar to multiply.
   * @param obj The IntegrandInter to multiply.
   * @return The resulting IntegrandInter after multiplication.
   */
  friend IntegrandInter operator*(const dcmplx& scalar, const IntegrandInter& obj) {
    IntegrandInter result = obj;
    multiplyVector(result.te,    scalar);
    multiplyVector(result.tm,    scalar);
    multiplyVector(result.tez1,  scalar);
    multiplyVector(result.tmz1,  scalar);
    multiplyVector(result.tez2,  scalar);
    multiplyVector(result.tmz2,  scalar);
    multiplyVector(result.tezz,  scalar);
    multiplyVector(result.tmzz,  scalar);
    multiplyVector(result.tes,   scalar);
    multiplyVector(result.tms,   scalar);
    multiplyVector(result.ter,   scalar);
    multiplyVector(result.tmr,   scalar);
    multiplyVector(result.terz1, scalar);
    multiplyVector(result.tmrz1, scalar);
    multiplyVector(result.terz2, scalar);
    multiplyVector(result.tmrz2, scalar);
    return result;
  }

  /**
   * @brief Function to create a flattened structure.
   * @param indVecTriplet The vector of index triplets.
   * @return The resulting flattened IntegrandInter.
   */
  IntegrandInter flatten(const std::vector<std::tuple<int, int, int>>& indVecTriplet) const {
    IntegrandInter result;

    // Process each member
    result.te    = flattenVector(te,    indVecTriplet);
    result.tm    = flattenVector(tm,    indVecTriplet);
    result.tez1  = flattenVector(tez1,  indVecTriplet);
    result.tmz1  = flattenVector(tmz1,  indVecTriplet);
    result.tez2  = flattenVector(tez2,  indVecTriplet);
    result.tmz2  = flattenVector(tmz2,  indVecTriplet);
    result.tezz  = flattenVector(tezz,  indVecTriplet);
    result.tmzz  = flattenVector(tmzz,  indVecTriplet);
    result.tes   = flattenVector(tes,   indVecTriplet);
    result.tms   = flattenVector(tms,   indVecTriplet);
    result.ter   = flattenVector(ter,   indVecTriplet);
    result.tmr   = flattenVector(tmr,   indVecTriplet);
    result.terz1 = flattenVector(terz1, indVecTriplet);
    result.tmrz1 = flattenVector(tmrz1, indVecTriplet);
    result.terz2 = flattenVector(terz2, indVecTriplet);
    result.tmrz2 = flattenVector(tmrz2, indVecTriplet);

    return result;
  }

  private:
  /**
   * @brief In-place element-wise subtraction of two 3D arrays.
   * @param[in,out] a Left-hand 3D array; updated with \c a-=b for overlapping extents.
   * @param[in] b     Right-hand 3D array to subtract from \p a.
   * @details Iterates the three nested dimensions and subtracts corresponding entries.
   *          If a size mismatch is found at any nesting level, that branch is skipped
   *          (or the whole operation returns early at the outermost level).
   * @note No allocations are performed; only existing elements in \p a are modified.
   */
  static void subtractVector(std::vector<std::vector<std::vector<dcmplx>>>& a,
                             const std::vector<std::vector<std::vector<dcmplx>>>& b) {
    if (a.size() != b.size()) return; // Ensure sizes match
    for (size_t i = 0; i < a.size(); ++i) {
      if (a[i].size() != b[i].size()) continue;
      for (size_t j = 0; j < a[i].size(); ++j) {
        if (a[i][j].size() != b[i][j].size()) continue;
        for (size_t k = 0; k < a[i][j].size(); ++k) {
          a[i][j][k] -= b[i][j][k];
        }
      }
    }
  }

  /**
   * @brief In-place scalar multiplication of a 3D array.
   * @param[in,out] vec   Target 3D array whose elements are scaled.
   * @param[in]     scalar Complex scaling factor applied to every element.
   * @details Traverses all dimensions and multiplies each entry of \p vec by \p scalar.
   * @note Operates in place; no copies are created.
   */
  static void multiplyVector(std::vector<std::vector<std::vector<dcmplx>>>& vec,
                             const dcmplx& scalar) {
    for (size_t i = 0; i < vec.size(); ++i) {
      for (size_t j = 0; j < vec[i].size(); ++j) {
        for (size_t k = 0; k < vec[i][j].size(); ++k) {
          vec[i][j][k] *= scalar;
        }
      }
    }
  }

  /**
   * @brief Extract selected elements from a 3D array into a compact 3D structure.
   * @param[in] matrix        Source 3D array.
   * @param[in] indVecTriplet List of index triplets \c (i,j,k) to extract.
   * @return A 3D array where each selected element is stored as a 1×1 "slice"
   *         (shape \c [[value]]), preserving 3D semantics for downstream code.
   * @details Each value is obtained via @ref getValue() (bounds-checked) and appended
   *          as a 1×1 block to the output.
   * @warning Out-of-bounds indices yield a zero value (see @ref getValue()).
   */
  std::vector<std::vector<std::vector<dcmplx>>> flattenVector(
      const std::vector<std::vector<std::vector<dcmplx>>>& matrix,
      const std::vector<std::tuple<int, int, int>>& indVecTriplet) const {
    std::vector<std::vector<std::vector<dcmplx>>> flattened;
    for (const auto& coord : indVecTriplet) {
      int i = std::get<0>(coord); // Extract first index
      int j = std::get<1>(coord); // Extract second index
      int k = std::get<2>(coord); // Extract third index
      flattened.push_back({ {getValue(matrix, i, j, k)} }); // Store as 1-column 3D vector
    }
    return flattened;
  }

  /**
   * @brief Bounds-checked access to a 3D array element.
   * @param[in] matrix Source 3D array.
   * @param[in] i      First (0-based) index.
   * @param[in] j      Second (0-based) index.
   * @param[in] k      Third (0-based) index.
   * @return The element at \c (i,j,k) if all indices are valid; otherwise \c dcmplx(0,0).
   * @note Never throws; returns a default complex for invalid indices.
   */
  dcmplx getValue(const std::vector<std::vector<std::vector<dcmplx>>>& matrix, int i, int j, int k) const {
    if (i >= 0 && i < static_cast<int>(matrix.size()) &&
        j >= 0 && j < static_cast<int>(matrix[i].size()) &&
        k >= 0 && k < static_cast<int>(matrix[i][j].size())) {
      return matrix[i][j][k];
    }
    return dcmplx(0, 0); // Default to zero if out of bounds
  }
};

/// \ingroup domain
/**
 * @class LayeredMediaUtils
 * @brief Utility class for calculations related to layered media.
 */
class LayeredMediaUtils {
  private:
    std::vector<dcmplx> k_L;              ///< Wavenumbers for layers.
    std::vector<dcmplx> epsilon;          ///< Permittivity values for layers.
    std::vector<dcmplx> mu;               ///< Permeability values for layers.
    std::vector<double> zValsInterfaces;  ///< Z-coordinates of layer interfaces.
    std::vector<double> thickness;        ///< Thickness of each layer.
    int midLayerIndex;                    ///< Layer index of the mesh middle point.

  public:
    /**
     * @brief Constructor initializing layered media parameters.
     * @param inLayeredK Vector of wavenumbers for each layer.
     * @param inEpsilonMedium Vector of relative permittivity values for each layer.
     * @param inMuMedium Vector of relative permeability values for each layer.
     * @param inZValsInterfaces Vector of Z-coordinates defining the layer interfaces.
     * @param inThickness Vector of thickness values for each layer.
     * @param inMidLayerIndex Layer index of the middle point in the mesh.
     */
    LayeredMediaUtils(const std::vector<dcmplx>& inLayeredK,
                      const std::vector<dcmplx>& inEpsilonMedium,
                      const std::vector<dcmplx>& inMuMedium,
                      const std::vector<double>& inZValsInterfaces,
                      const std::vector<double>& inThickness,
                      int inMidLayerIndex = 0);
 
    /** @brief Destructor. */
    ~LayeredMediaUtils() = default;

    /**
     * @brief Returns the layer index for a given spatial position.
     * @param pos The spatial position as an rvec (3D vector).
     * @return The layer index (0 to N-1), where N is the number of layers.
     */
    int GetLayerIndex(rvec pos) const;

    /**
     * @brief Computes Fresnel reflection and transmission coefficients at an interface.
     * @param layerIndexI Index of the first layer.
     * @param layerIndexJ Index of the second layer.
     * @param krho Transverse wavenumber.
     * @param waves Type of wave polarization ("TE" or "TM").
     * @return A vector of complex values representing the Fresnel coefficients.
     */
    std::vector<dcmplx> FresnelCoeff(int layerIndexI, int layerIndexJ, 
                                     dcmplx krho, const std::string& waves) const;                                  
    
    /**
     * @brief Computes the interface transfer matrix for TE or TM waves.
     * @param layerIndexI Index of the first layer.
     * @param layerIndexJ Index of the second layer.
     * @param krho Transverse wavenumber.
     * @param waves Type of wave polarization ("TE" or "TM").
     * @return A 2D vector of complex values representing the interface transfer matrix.
     */
    std::vector<std::vector<dcmplx>> InterfaceTransferMatrix(
                                      int layerIndexI, int layerIndexJ,
                                      dcmplx krho, const std::string& waves) const;
    
    /**
     * @brief Computes the propagation matrix for a given layer.
     * @param layerIndex Index of the layer for which the propagation matrix is computed.
     * @param krho Transverse wavenumber.
     * @return A 2x2 matrix representing the propagation matrix of the layer.
     * @throws std::out_of_range if the layer index is invalid.
     */
    std::vector<std::vector<dcmplx>> propagationMatrix(int layerIndex, dcmplx krho) const;
    
    /**
     * @brief Computes the total transfer matrix for the entire layered structure.
     * @param krho Transverse wavenumber.
     * @param waves Type of wave polarization ("TE" or "TM").
     * @return A 2x2 matrix representing the total transfer matrix of the layered medium.
     * @throws std::invalid_argument if an invalid wave type is provided.
     */
    std::vector<std::vector<dcmplx>> totalTransferMatrix(dcmplx krho, const std::string& waves) const;

    /**
     * @brief Computes reflection and transmission coefficients for the entire layered structure.
     * @param krho Transverse wavenumber.
     * @param waves Type of wave polarization ("TE" or "TM").
     * @param direction Direction of wave propagation ("up" or "down").
     * @return A vector containing the total Fresnel reflection and transmission coefficients.
     * @throws std::invalid_argument if an invalid wave type or direction is provided.
     */
    std::vector<dcmplx> totalFresnelCoeffs(dcmplx krho, const std::string& waves, 
                                                        const std::string& direction) const;
    
    /**
     * @brief Computes secondary wave coefficients for a given observation and source layer.
     * @param krho Transverse wavenumber.
     * @param observationLayer Index of the observation layer.
     * @param sourceLayer Index of the source layer.
     * @param waves Type of wave polarization ("TE" or "TM").
     * @return A 2D Blitz array of complex values representing the secondary wave coefficients.
     */
    blitz::Array<dcmplx, 2> Secondary(dcmplx krho, int observationLayer, 
                                      int sourceLayer, const std::string& waves);
    
    /**
     * @brief Computes electromagnetic fields at observation points given a source position.
     * @param krho Transverse wavenumber.
     * @param observationPointsZ Vector of Z-coordinates for observation points.
     * @param sourcePointZ Z-coordinate of the source point.
     * @param waves Type of wave polarization ("TE" or "TM").
     * @return A vector of complex values representing the computed field values.
     */
    std::vector<dcmplx> Fields(dcmplx krho, std::vector<double> observationPointsZ, 
                               double sourcePointZ, const std::string& waves);

    /**
     * @brief Performs matrix multiplication for 2D complex matrices.
     * @param x The first matrix (left operand).
     * @param y The second matrix (right operand).
     * @return The resulting matrix after multiplication.
     * @throws std::invalid_argument if the dimensions of the matrices are incompatible.
     */
    std::vector<std::vector<dcmplx>> matrixMultiply2D(const std::vector<std::vector<dcmplx>>& x,
                                                      const std::vector<std::vector<dcmplx>>& y) const;

    /**
     * @brief Performs matrix multiplication for 3D complex matrices.
     * @param x The first 3D matrix (left operand).
     * @param y The second 3D matrix (right operand).
     * @return The resulting 3D matrix after element-wise multiplication.
     */
    std::vector<std::vector<std::vector<dcmplx>>> matrixMultiply3D(
                                                      const std::vector<std::vector<std::vector<dcmplx>>>& x,
                                                      const std::vector<std::vector<std::vector<dcmplx>>>& y) const;

    /**
     * @brief Performs matrix multiplication for 2D/3D complex matrices.
     * @param x The first matrix (2D or 3D).
     * @param y The second matrix (2D or 3D).
     */
    void matrixMultiplication(const std::vector<std::vector<std::vector<dcmplx>>>& x,
                              const std::vector<std::vector<std::vector<dcmplx>>>& y) const;
    
    /**
     * @brief Transforms Cartesian coordinates to stratified medium coordinates
     * for the calculation of the quasistatic contribution of Green's functions.
     * @param r1 Observation point in Cartesian coordinates (x, y, z).
     * @param r2 Source point in Cartesian coordinates (x, y, z).
     * @param type Type of interaction: "intraTB" (intra-layer for top/bottom layer), 
     * "intraIn" (intra-layer for interior layer), or "inter" (inter-layer).
     * @return A layeredCoords struct containing the converted coordinates.
     */
    layeredCoords cartesianToLayeredQS(const rvec &r1, const rvec &r2, 
                                       const std::string& type) const;
    
    /**
     * @brief Transforms Cartesian coordinates to stratified medium coordinates
     * for the smooth Green's functions in the tabulation.
     * @param r1 Observation point in Cartesian coordinates (x, y, z).
     * @param r2 Source point in Cartesian coordinates (x, y, z).
     * @param type Type of interaction: "intraTB" (intra-layer for top/bottom layer),
     * "intraIn" (intra-layer for interior layer), or "inter" (inter-layer).
     * @return A layeredCoords struct containing the converted coordinates.
     */
    layeredCoords cartesianToLayeredTab(const rvec &r1, const rvec &r2, 
                                        const std::string& type) const;

    /**
     * @brief Computes the sign function for a given value.
     * @param val The value to compute the sign for.
     * @return 1 if val > 0, -1 if val < 0, and 0 if val == 0.
     */
    template <typename T> int sgn(T val) {
      return (T(0) < val) - (val < T(0));
    }

    /**
     * @brief Getter of k_L.
     * @param index The index of the k_L value to retrieve.
     * @return The k_L value at the specified index.
     */
    dcmplx getkL(int index) {
      return k_L[index];
    }

    /**
     * @brief Getter of zValInterface.
     * @param index The index of the zValInterface value to retrieve.
     * @return The zValInterface value at the specified index.
     */
    double getZValInterface(int index) {
      if (index < 0 || index >= zValsInterfaces.size()) {
        throw std::out_of_range("Index out of range in getZValInterface");
      }
      return zValsInterfaces[index];
    }

    /**
     * @brief Getter of epsilon.
     * @param index The index of the epsilon value to retrieve.
     * @return The epsilon value at the specified index.
     */
    dcmplx getEps(int index) {
      if (index < 0 || index >= epsilon.size()) {
        throw std::out_of_range("Index out of range in getEpsilon");
      }
      return epsilon[index];
    }

    /**
     * @brief Getter of mu.
     * @param index The index of the mu value to retrieve.
     * @return The mu value at the specified index.
     */
    dcmplx getMu(int index) {
      if (index < 0 || index >= mu.size()) {
        throw std::out_of_range("Index out of range in getMu");
      }
      return mu[index];
    }

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
};

#endif
