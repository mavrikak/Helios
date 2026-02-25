/**
 * @file mathFunctions.h
 * @brief Header file of the class that provides auxiliary mathematical functionalities
 * for the SIE solver. These include:
 * - Lookup table management classes for complex-valued functions
 * - Special functions (erf/erfc, exponential integral)
 * - Factorials and dyadic algebra
 * - Symmetry operations (translations)
 * - Interpolation templates (linear, bilinear, trilinear)
 * - Table generators for error-function lookups (text and binary)
 */

#ifndef MATHFUNCTIONS_H
#define MATHFUNCTIONS_H

#include <math.h>
#include <complex>
#include <vector>
#include "globals.h"

/// \ingroup greenFunction
/**
 * @class LookupTable
 * @brief Represents a 2D lookup table for complex functions (text format).
 *
 * Stores grid of sampled argument (arg) and function values (fun).
 */
class LookupTable {
 public:
  std::vector<std::vector<dcmplx>> arg; ///< Argument grid (Re,Im)
  std::vector<std::vector<dcmplx>> fun; ///< Function values
  int sizeRe;                           ///< Number of samples in real direction
  int sizeIm;                           ///< Number of samples in imaginary direction
  double maxRe;                         ///< Maximum ranges for real/imag axes
  double maxIm;                         ///< Maximum ranges for real/imag axes
  double incRe;                         ///< Sampling increments
  double incIm;                         ///< Sampling increments

  // ---------------------------------------------------------------------------
  // Construction & lifetime
  // ---------------------------------------------------------------------------
  /**
   * @brief Default constructor. Initializes table parameters to zero.
   */
  LookupTable();

  /**
   * @brief Construct a LookupTable from a text file.
   * @param fileName Path to the lookup table text file.
   * @details Reads sampled argument and function values into memory.
   */
  LookupTable(std::string fileName);
};

/// \ingroup greenFunction
/**
 * @class LookupTableBin
 * @brief Represents a 2D lookup table for complex functions (binary format).
 */
class LookupTableBin {
 public:
  std::vector<std::vector<dcmplx>> arg; ///< Argument grid (Re,Im)
  std::vector<std::vector<dcmplx>> fun; ///< Function values
  int sizeRe;                           ///< Number of samples in real direction
  int sizeIm;                           ///< Number of samples in imaginary direction
  double maxRe;                         ///< Maximum ranges for real/imag axes
  double maxIm;                         ///< Maximum ranges for real/imag axes
  double incRe;                         ///< Sampling increments
  double incIm;                         ///< Sampling increments

  // ---------------------------------------------------------------------------
  // Construction & lifetime
  // ---------------------------------------------------------------------------
  /**
   * @brief Default constructor. Initializes table parameters to zero.
   */
  LookupTableBin();

  /**
   * @brief Construct a LookupTableBin from a binary file.
   * @param fileName Path to the lookup table binary file.
   * @details Reads sampled argument and function values into memory.
   */
  LookupTableBin(std::string fileName);
};

// -----------------------------------------------------------------------------
// Symmetry operations
// -----------------------------------------------------------------------------
/// \ingroup greenFunction
/**
 * @class SymmetryOperation
 * @brief Abstract base class for spatial symmetry operations.
 */
class SymmetryOperation {
 public:
  /** @brief Default destructor */
  virtual ~SymmetryOperation(){};

  /**
   * @brief Apply a translation to a given vector.
   * @param inputVector Input position vector.
   * @return Translated vector (input minus translation vector).
   */
  virtual rvec Apply(rvec inputVector) = 0;
};

/// \ingroup greenFunction
/**
 * @class Translation
 * @brief Concrete translation symmetry operation.
 */
class Translation : public SymmetryOperation {
 private:
  rvec t;   ///< Translation vector

 public:
  /**
   * @brief Constructor.
   * @param translationVector Translation vector to apply to coordinates.
   */
  Translation(rvec translationVector);

  /**
   * @brief Apply the translation to a given vector.
   * @param inputVector Input position vector.
   * @return Translated vector (input minus translation vector).
   */
  virtual rvec Apply(rvec inputVector);
};

// -----------------------------------------------------------------------------
// Special functions
// -----------------------------------------------------------------------------
/** 
 * @brief Error function for complex arguments.
 * @param arg Complex argument of \p cerf().
 * @return Complex result.
 */
dcmplx cerf(dcmplx arg);

/** 
 * @brief Complementary error function for complex arguments.
 * @param arg Complex argument of \p cerfc().
 * @return Complex result.
 */
dcmplx cerfc(dcmplx arg);

/** 
 * @brief Exponential integral function \p E1/Ei(x).
 * @param x Real argument of \p cerfc().
 * @return Real result.
 */
double Ei(const double &x);

/** 
 * @brief Exponential integral of @f$ p^{th} @f$ order.
 * @param x Real argument of integral.
 * @param p Order of integral.
 * @return Complex result.
 */
dcmplx Ep(const double &x, int p);

/** 
 * @brief Factorial.
 * @param num Integer argument of factorial.
 * @return Integer result.
 */
int factorial(int num);

/**
 * @brief Complex 3x3 dyadic and dot products.
 * @param d Complex dyadic.
 * @return Transposed complex dyadic.
 */
cdyad transpose(cdyad d);

/**
 * @brief Evaluate a Green's function at given source and observation points.
 * @tparam Type Return type of the Green's function (e.g., complex or scalar).
 * @param r Observation point.
 * @param rp Source point.
 * @param greenFunction Pointer to the Green's function to evaluate.
 * @details Prints the evaluated Green's function value to standard output.
 */
template <class Type>
void FillGreenTable(rvec r, rvec rp, Type (*greenFunction)(rvec, rvec)) {
  std::cout << greenFunction(r, rp) << std::endl;
}

// -----------------------------------------------------------------------------
// Interpolation templates
// -----------------------------------------------------------------------------
/**
 * @brief Perform linear interpolation between two samples.
 * @tparam Type Data type of interpolated values.
 * @param x0 First sample position.
 * @param i0 First sample value.
 * @param x1 Second sample position.
 * @param i1 Second sample value.
 * @param xp Query position between x0 and x1.
 * @return Interpolated value at xp.
 */
template <class Type>
Type InterpLinear(double x0, Type i0, double x1, Type i1, float xp) {
  return (x1 - xp) / (x1 - x0) * i0 + (xp - x0) / (x1 - x0) * i1;
}

/**
 * @brief Perform bilinear interpolation on a rectangular surface.
 * @tparam Type Data type of interpolated values.
 * @param x0,y0,i0 Coordinates and value of first vertex.
 * @param x1,i1 Coordinates and value of second vertex.
 * @param x2,i2 Coordinates and value of third vertex.
 * @param x3,y3,i3 Coordinates and value of fourth vertex.
 * @param xp,yp Query coordinates inside the quadrilateral.
 * @return Interpolated value at (xp, yp).
 */
template <class Type>
Type InterpBilinear(double x0, double y0, Type i0,
                    double x1, double /*y1*/, Type i1,
                    double x2, double /*y2*/, Type i2, 
                    double x3, double y3, Type i3,
                    double xp, double yp) {
  // clang-format on
  Type i01 = InterpLinear<Type>(x0, i0, x1, i1, xp);
  Type i32 = InterpLinear<Type>(x3, i3, x2, i2, xp);

  return InterpLinear<Type>(y0, i01, y3, i32, yp);
}

/**
 * @brief Perform trilinear interpolation inside a cube.
 * @tparam Type Data type of interpolated values.
 * @param x0, y0, z0, i0 Coordinates and value of the 1st cube corner.
 * @param x1, y1, i1 Coordinates and value of the 2nd cube corner.
 * @param x2, y2, i2 Coordinates and value of the 3rd cube corner.
 * @param x3, y3, i3 Coordinates and value of the 4th cube corner.
 * @param x4, y4, z4, i4 Coordinates and value of the 5th cube corner.
 * @param x5, y5, i5 Coordinates and value of the 6th cube corner.
 * @param x6, y6, i6 Coordinates and value of the 7th cube corner.
 * @param x7, y7, i7 Coordinates and value of the 8th cube corner.
 * @param xp,yp,zp Query coordinates inside the cube.
 * @return Interpolated value at (xp, yp, zp).
 */
template <class Type>
Type InterpTrilinear(double x0, double y0, double z0, Type i0,
                     double x1, double y1, double /*z1*/, Type i1,
                     double x2, double y2, double /*z2*/, Type i2,
                     double x3, double y3, double /*z3*/, Type i3,
                     double x4, double y4, double z4, Type i4, 
                     double x5, double y5, double /*z5*/, Type i5,
                     double x6, double y6, double /*z6*/, Type i6,
                     double x7, double y7, double /*z7*/, Type i7,
                     double xp, double yp, double zp) {
  // clang-format on
  Type i0123 = InterpBilinear<Type>(x0, y0, i0, x1, y1, i1, x2, y2, i2, x3, y3,
                                    i3, xp, yp);

  Type i4567 = InterpBilinear<Type>(x4, y4, i4, x5, y5, i5, x6, y6, i6, x7, y7,
                                    i7, xp, yp);

  return InterpLinear<Type>(z0, i0123, z4, i4567, zp);
}

// -----------------------------------------------------------------------------
// Table generation helpers
// -----------------------------------------------------------------------------
/**
 * @brief Generate and save a text lookup table for \p cerfc over a complex grid.
 * @param filename Output text file name.
 * @param maxRe Maximum real value of the grid.
 * @param maxIm Maximum imaginary value of the grid.
 * @param incRe Increment along the real axis.
 * @param incIm Increment along the imaginary axis.
 * @return 0 on success, nonzero on failure.
 */
int ErfcLookupTable(std::string filename, double maxRe, double maxIm,
                    double incRe, double incIm);

/**
 * @brief Generate and save a binary lookup table for \p cerfc over a complex grid.
 * @param filename Output binary file name.
 * @param maxRe Maximum real value of the grid.
 * @param maxIm Maximum imaginary value of the grid.
 * @param incRe Increment along the real axis.
 * @param incIm Increment along the imaginary axis.
 * @return 0 on success, nonzero on failure.
 */
int ErfcLookupTableBin(std::string filename, double maxRe, double maxIm,
                       double incRe, double incIm);

#endif
