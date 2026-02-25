/**
 * @file SommerfeldIntegrator.h
 * @brief Declaration of the SommerfeldIntegrator class, which evaluates
 *        Sommerfeld-type integrals for layered media Green's functions.
 *
 * The class provides:
 *  - Path selection and data packing for semi-ellipse, real-inf, and imag-inf contours.
 *  - Intra-/inter-layer integrand assembly (TE/TM, field/derivative/surface terms).
 *  - Numerical integration with singularity subtraction (Chew, W.C., Tong, M.S., 
 *    Hu, B. (2009), Chapter 6).
 *  - Interpolation helpers to compact results for later evaluation.
 *
 * Types like #Data2D/#Data3D, #IntegrandIntra/#IntegrandInter, and 
 * #CoeffsIntra/#CoeffsInter are defined in headers included transitively 
 * (e.g., LayeredMediaUtils.h).
 */

#ifndef SOMMERFELD_INTEGRATOR_H
#define SOMMERFELD_INTEGRATOR_H

#include <map>
#include "DomainHom3D.h"
#include "SurfaceMesh.h"
#include "LayeredMediaUtils.h"
#include <stdexcept>         // For exception handling
#include <functional>        // For std::function

/// \ingroup greenFunction
/**
 * @class SommerfeldIntegrator
 * @brief Numerical engine for Sommerfeld integrals in stratified media.
 *
 * Construct with material stacks and numerical knobs (semi-ellipse radii, ratio, cutoff),
 * then use Evaluate() to get intra-/inter-layer Green coefficients packed for interpolation.
 */
class SommerfeldIntegrator {
  private:
    // ---- Material stack / geometry ----
    std::vector<dcmplx> k_L;               ///< Wavevector of the layered medium
    int materialIndex;                     ///< Observation  layer index (and source for intra-layer)
    int sourceIndex;                       ///< Source layer index
    std::vector<dcmplx> epsilon;           ///< Permittivity of the layered medium
    std::vector<dcmplx> mu;                ///< Permeability of the layered medium
    std::vector<double> zValsInterfaces;   ///< Interface z-values
    std::vector<double> thickness;         ///< Layer thicknesses

    // ---- Contour / subtraction parameters ----
    double semiReal;                       ///< Scaling factor for real axis of semiellipse
    double semiImag;                       ///< Scaling factor for imaginary axis of semiellipse
    double ratio;                          ///< z:r ratio
    std::map<std::string, double> options; ///< Options for ODE integration
    double cutoff;                         ///< Cutoff parameter

    // ---- Utilities ----
    LayeredMediaUtils* layeredUtils;       ///< Layered media utility functions

  public:
    /**
     * @brief Constructs a SommerfeldIntegrator object for layered media calculations.
     *
     * This constructor initializes the Sommerfeld integrator with the given parameters,
     * which define the properties of the layered medium and numerical settings for integration.
     *
     * @param k_L A vector of complex wave numbers for the layered medium.
     * @param inMaterialIndex An integer indicating the observation/source intra-layer index for calculations.
     * @param inEpsilonMedium A vector of complex permittivity values for the medium layers.
     * @param inMuMedium A vector of complex permeability values for the medium layers.
     * @param inZValsInterfaces A vector of interface positions along the z-axis.
     * @param inThickness A vector of layer thicknesses.
     * @param inSemiReal Scaling factor for real axis of semiellipse.
     * @param inSemiImag Scaling factor for imaginary axis of semiellipse.
     * @param inRatio A double defining the integration ratio parameter.
     * @param inCutoff A double specifying the cutoff value for the integration.
     */
    SommerfeldIntegrator(const std::vector<dcmplx> k_L, int inMaterialIndex,
                         const std::vector<dcmplx>& inEpsilonMedium,
                         const std::vector<dcmplx>& inMuMedium,
                         const std::vector<double>& inZValsInterfaces,
                         const std::vector<double>& inThickness,
                         double inSemiReal = 4.0, double inSemiImag = 0.05,
                         double inRatio = 1.0, double inCutoff = 0.1);
    
    /**
     * @brief Constructs a SommerfeldIntegrator object for layered media calculations.
     *
     * This constructor initializes the Sommerfeld integrator with the given parameters,
     * which define the properties of the layered medium and numerical settings for integration.
     *
     * @param k_L A vector of complex wave numbers for the layered medium.
     * @param inMaterialIndex An integer indicating the observation layer index for calculations.
     * @param inSourceIndex An integer indicating the source layer index for calculations.
     * @param inEpsilonMedium A vector of complex permittivity values for the medium layers.
     * @param inMuMedium A vector of complex permeability values for the medium layers.
     * @param inZValsInterfaces A vector of interface positions along the z-axis.
     * @param inThickness A vector of layer thicknesses.
     * @param inSemiReal Scaling factor for real axis of semiellipse.
     * @param inSemiImag Scaling factor for imaginary axis of semiellipse.
     * @param inRatio A double defining the integration ratio parameter.
     * @param inCutoff A double specifying the cutoff value for the integration.
     */
    SommerfeldIntegrator(const std::vector<dcmplx> k_L, 
                         int inMaterialIndex, int inSourceIndex, 
                         const std::vector<dcmplx>& inEpsilonMedium,
                         const std::vector<dcmplx>& inMuMedium,
                         const std::vector<double>& inZValsInterfaces,
                         const std::vector<double>& inThickness,
                         double inSemiReal = 4.0, double inSemiImag = 0.05,
                         double inRatio = 0.8, double inCutoff = 0.1);
 
    /** @brief Destructor: releases LayeredMediaUtils. */
    ~SommerfeldIntegrator();
    
    // -------------------- Path selection helpers --------------------
    /**
     * @brief Selects and processes 2D data based on the given input parameters.
     *
     * This function categorizes and stores 2D data depending on the provided path mode ('semi', 'real', or 'imag').
     * It organizes radial and axial coordinates into matrix or vector formats as needed.
     *
     * @param r A 2D vector containing the radial coordinate values.
     * @param z1 A 2D vector containing the first axial coordinate values.
     * @param z2 A 2D vector containing the second axial coordinate values (optional).
     * @param path A string specifying the mode of data selection ('semi', 'real', or 'imag').
     *
     * @return A Data2D structure containing processed coordinate values and indexing information.
     *
     * @throws std::invalid_argument if an invalid path is specified.
     */
    Data2D selectData2D(const std::vector<std::vector<double>>& r,
                        const std::vector<std::vector<double>>& z1,
                        const std::vector<std::vector<double>>& z2,
                        const std::string& path);
    
    /**
     * @brief Selects and processes 3D data based on the given input parameters.
     *
     * This function categorizes and stores 3D data depending on the provided path mode ('semi', 'real', or 'imag').
     * It organizes radial and axial coordinates into matrix or vector formats as needed.
     *
     * @param r A 3D vector containing the radial coordinate values.
     * @param z1 A 3D vector containing the first axial coordinate values.
     * @param z2 A 3D vector containing the second axial coordinate values.
     * @param path A string specifying the mode of data selection ('semi', 'real', or 'imag').
     *
     * @return A Data3D structure containing processed coordinate values and indexing information.
     *
     * @throws std::invalid_argument if an invalid path is specified.
     */
    Data3D selectData3D(const std::vector<std::vector<std::vector<double>>>& r,
                        const std::vector<std::vector<std::vector<double>>>& z1,
                        const std::vector<std::vector<std::vector<double>>>& z2,
                        const std::string& path);

    /**
     * @brief Computes distances from interfaces for inter-layer secondary wave calculations.
     *
     * This function calculates the distances of given points from an interface of a layer, 
     * adjusting the calculations based on the specified index and direction.
     *
     * @param z A 3D vector representing the z-coordinates of points.
     * @param index The index of the layer of interest.
     * @param direction A double indicating the direction of propagation.
     *
     * @return A 3D vector containing the computed distances from the interface.
     */
    std::vector<std::vector<std::vector<double>>> interDist(
                                            std::vector<std::vector<std::vector<double>>> z, 
                                            int index, double direction);
    
    /**
     * @brief Computes secondary wave contributions in a layered medium.
     *
     * This function calculates secondary wave fields based on given wavenumber, 
     * radial and axial coordinates, and the type of waves involved.
     *
     * @param krho A complex wavenumber in the radial direction.
     * @param r A 3D vector containing the radial coordinate values.
     * @param z1 A 3D vector containing the first axial coordinate values.
     * @param z2 A 3D vector containing the second axial coordinate values.
     * @param waves A string specifying the type of waves considered (e.g., TE, TM).
     *
     * @return A 4D vector containing:
     *         - f: The computed field values.
     *         - fZ1: The first Z1-derivative of the field.
     *         - fZ2: The first Z2-derivative of the field.
     *         - fZZ: The mixed Z1-Z2-derivative of the field.
     */
    std::vector<std::vector<std::vector<std::vector<dcmplx>>>> 
    Secondary(dcmplx krho, std::vector<std::vector<std::vector<double>>> r,
              std::vector<std::vector<std::vector<double>>> z1,
              std::vector<std::vector<std::vector<double>>> z2, const std::string& waves);
    
    // -------------------- Integrands --------------------
    /**
     * @brief Computes the intra-layer Green's function integrand.
     *
     * This function evaluates the intra-layer integrand for a given wavenumber, using radial and axial coordinates.
     * It computes Green function elements, derivatives, and surface integrals based on mode selection and operations.
     *
     * @param rq_te Quasistatic TE reflection coefficients.
     * @param rq_tm Quasistatic TM reflection coefficients.
     * @param inData2D Data2D structure containing spatial coordinate information.
     * @param kr Complex radial wavenumber.
     * @param kz Complex axial wavenumber.
     * @param mode String specifying the mode of operation ('bessel' or 'hankel').
     * @param op String specifying additional operations ('plus', 'minus', or empty for default behavior).
     *
     * @return An IntegrandIntra structure containing computed Green function elements and derivatives.
     *
     * @throws std::runtime_error if input data dimensions mismatch or invalid mode/op is provided.
     */
    IntegrandIntra IntraIntegrand(const std::vector<dcmplx>& rq_te, 
                                  const std::vector<dcmplx>& rq_tm,
                                  const Data2D& inData2D, dcmplx kr, dcmplx kz, 
                                  const std::string& mode, const std::string& op = "");
    
    /**
     * @brief Computes the inter-layer Green's function integrand.
     *
     * This function evaluates the inter-layer integrand for a given wavenumber, using radial and axial coordinates.
     * It computes Green function elements, derivatives, and surface integrals based on mode selection and operations.
     *
     * @param tq_te Quasistatic TE transmission coefficient.
     * @param tq_tm Quasistatic TM transmission coefficient.
     * @param inData3D Data3D structure containing spatial coordinate information.
     * @param kr Complex radial wavenumber.
     * @param mode String specifying the mode of operation ('bessel' or 'hankel').
     *
     * @return An IntegrandInter structure containing computed Green function elements and derivatives.
     *
     * @throws std::runtime_error if input data dimensions mismatch or invalid mode/op is provided.
     */
    IntegrandInter InterIntegrand(const dcmplx& tq_te, const dcmplx& tq_tm,
                                  const Data3D& inData3D, dcmplx kr, 
                                  const std::string& mode);
    
    // -------------------- Integration and evaluation --------------------
    /**
     * @brief Computes the intra-layer integration function.
     *
     * This function evaluates the integral of the intra-layer Green's function over a specified path.
     * The integration can be performed along a semi-ellipse, real infinity, or imaginary infinity.
     *
     * @param rq_te Quasistatic TE reflection coefficients.
     * @param rq_tm Quasistatic TM reflection coefficients.
     * @param inData A Data2D structure containing spatial coordinate information.
     * @param x Double representing the integration parameter.
     * @param k_max Parameter for integration path specification.
     * @param intZeroKrho An IntegrandIntra structure containing precomputed singularity corrections.
     * @param path A string specifying the integration path ('semi', 'real', or 'imag').
     * @param op A string specifying additional operations ('plus', 'minus', or empty for default behavior).
     *
     * @return An IntegrandIntra structure containing the integrated Green function elements.
     *
     * @throws std::invalid_argument if an invalid path type is specified.
     */
    IntegrandIntra integrationFIntra(const std::vector<dcmplx>& rq_te, 
                                     const std::vector<dcmplx>& rq_tm, 
                                     const Data2D& inData, double x, double k_max,
                                     const IntegrandIntra& intZeroKrho, 
                                     const std::string& path, const std::string& op = "");
    
    /**
     * @brief Computes the inter-layer integration function.
     *
     * This function evaluates the integral of the inter-layer Green's function over a specified path.
     * The integration can be performed along a semi-ellipse, real infinity, or imaginary infinity.
     *
     * @param tq_te Quasistatic TE transmission coefficient.
     * @param tq_tm Quasistatic TM transmission coefficient.
     * @param inData A Data3D structure containing spatial coordinate information.
     * @param x Double representing the integration parameter.
     * @param k_max Parameter for integration path specification.
     * @param intZeroKrho An IntegrandInter structure containing precomputed singularity corrections.
     * @param path A string specifying the integration path ('semi', 'real', or 'imag').
     *
     * @return An IntegrandInter structure containing the integrated Green function elements.
     *
     * @throws std::invalid_argument if an invalid path type is specified.
     */
    IntegrandInter integrationFInter(const dcmplx& tq_te, const dcmplx& tq_tm, 
                                     const Data3D& inData, double x, double k_max,
                                     const IntegrandInter& intZeroKrho, 
                                     const std::string& path);
    
    /**
     * @brief Computes intra-layer coefficients for Green's function interpolation.
     *
     * @param rq_te Quasistatic TE reflection coefficients.
     * @param rq_tm Quasistatic TM reflection coefficients.
     * @param r A 2D vector representing radial coordinate values.
     * @param z A 2D vector representing axial coordinate values.
     * @param singular A boolean indicating whether singularity subtraction is required (always true).
     * @param op A string specifying additional operations ('plus', 'minus', or empty for default behavior).
     *
     * @return A CoeffsIntra structure containing Green's functions interpolation coefficients.
     *
     * @throws std::invalid_argument if an invalid singular value is provided.
     */
    CoeffsIntra Evaluate(const std::vector<dcmplx>& rq_te, const std::vector<dcmplx>& rq_tm,
                         const std::vector<std::vector<double>>& r, 
                         const std::vector<std::vector<double>>& z, 
                         bool singular, const std::string& op = "");
    
    /**
     * @brief Computes inter-layer coefficients for Green's function interpolation.
     *
     * @param tq_te Quasistatic TE transmission coefficient.
     * @param tq_tm Quasistatic TM transmission coefficient.
     * @param r A 3D vector representing radial coordinate values.
     * @param z1 A 3D vector representing axial coordinate values for the first layer.
     * @param z2 A 3D vector representing axial coordinate values for the second layer.
     * @param singular A boolean indicating whether singularity subtraction is required (always true).
     * @param op A string specifying additional operations (always empty).
     *
     * @return A CoeffsInter structure containing Green's functions interpolation coefficients.
     *
     * @throws std::invalid_argument if an invalid singular value is provided.
     */
    CoeffsInter Evaluate(const dcmplx& tq_te, const dcmplx& tq_tm,
                         const std::vector<std::vector<std::vector<double>>>& r, 
                         const std::vector<std::vector<std::vector<double>>>& z1, 
                         const std::vector<std::vector<std::vector<double>>>& z2, 
                         bool singular, const std::string& op = "");

    // -------------------- Debugging / testing --------------------
    /** 
     * @brief Dump a 2D complex matrix to text (i j Re Im) for debugging. 
     * @param matrix The 2D matrix to dump.
     * @param filename The name of the file to write to.
     */
    void writeMatrixToFile(const std::vector<std::vector<dcmplx>>& matrix, 
                           const std::string& filename);

    /** 
     * @brief Dump a 3D complex matrix to text (i j k Re Im) for debugging. 
     * @param matrix The 3D matrix to dump.
     * @param filename The name of the file to write to.
     */
    void write3DMatrixToFile(const std::vector<std::vector<std::vector<dcmplx>>>& matrix, 
                             const std::string& filename);
    
    // -------------------- Fresnel coefficients --------------------
    /**
     * @brief Computes the Fresnel reflection and transmission coefficients between two layers.
     *
     * This function calculates the Fresnel reflection and transmission coefficients
     * for TE (transverse electric) or TM (transverse magnetic) waves at the interface
     * between two material layers.
     *
     * @param layerIndexI Index of the first layer.
     * @param layerIndexJ Index of the second layer.
     * @param krho The transverse wave vector component.
     * @param waves A string specifying the wave type ('TE' or 'TM').
     *
     * @return A vector containing:
     *         - The Fresnel reflection coefficient.
     *         - The Fresnel transmission coefficient.
     *
     * @throws std::out_of_range if the layer indices are invalid.
     * @throws std::invalid_argument if the wave type is not 'TE' or 'TM'.
     */
    std::vector<dcmplx> FresnelCoeff(int layerIndexI, int layerIndexJ, 
                                     dcmplx krho, const std::string& waves) const;
};

#endif // SOMMERFELD_INTEGRATOR_H
