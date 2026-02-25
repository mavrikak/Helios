/**
 * @file PlaneWave.h
 * @brief Defines the PlaneWave class: representation of a plane‐wave incident field. 
 * Supports both homogeneous and layered media cases. Provides field evaluation,
 * integration with basis functions, and data export.
 */

#ifndef PLANEWAVE_H
#define PLANEWAVE_H

#include "GaussQuad.h"
#include "IncidentField.h"
#include "globals.h"
#include "LayeredMediaUtils.h"

/// \ingroup incidentField
/**
 * @class PlaneWave
 * @brief Models a plane wave as an incident field in homogeneous or layered media.
 */
class PlaneWave : public IncidentField {
 private:
  cvec k;                           ///< Wavevector
  dcmplx krho;                      ///< Transverse wavenumber component
  cvec p;                           ///< Polarization vector
  double z;                         ///< Impedance (homogeneous medium)
  dcmplx zL;                        ///< Impedance in layered medium
  rvec dir;                         ///< Propagation direction (unit vector)
  double incLayerZ;                 ///< z of first interface in incidence layer
  int incLayerIndex;                ///< Layer index of incidence
  GaussQuad<PlaneWave> gauss;       ///< Gaussian quadrature integrator
  LayeredMediaUtils* layeredUtils;  ///< Utilities for layered case

 public:
  /**
   * @brief Construct a homogeneous plane wave.
   * @param wvl Incident wavelength.
   * @param eps Permittivity of the homogeneous medium.
   * @param dir Propagation direction (3D vector).
   * @param po Polarization vector.
   */
  PlaneWave(dcmplx wvl, double eps, rvec dir, cvec po);
    
  /**
   * @brief Construct a plane wave in a layered medium.
   * @param k_L Wave numbers per layer.
   * @param epsL Permittivity values per layer.
   * @param muL Permeability values per layer.
   * @param zValsInterfaces z-coordinates of interfaces.
   * @param thickness Thicknesses of layers.
   * @param propDir Propagation direction.
   * @param po Polarization vector.
   */
  PlaneWave(std::vector<dcmplx> k_L, std::vector<dcmplx> epsL,
            std::vector<dcmplx> muL, std::vector<double> zValsInterfaces,
            std::vector<double> thickness, rvec propDir, cvec po);
  
  /** @brief Destructor. Cleans up LayeredMediaUtils if allocated. */
  ~PlaneWave() { if (layeredUtils) { delete layeredUtils; } }
  
  /**
   * @brief Evaluate the electric field at a given position (homogeneous).
   * @param r Observation point.
   * @return Electric field vector at the point.
   */
  cvec EvaluateE(rvec r);

  /**
   * @brief Evaluate the electric and magnetic fields at a given position in layered media.
   * @param r Observation point.
   * @return Vector containing electric and magnetic field complex vectors at the point.
   */
  std::vector<cvec> EvaluateLayeredEandH(rvec r);

  /**
   * @brief Evaluate the magnetic field at a given position (homogeneous).
   * @param r Observation point.
   * @return Magnetic field vector at the point.
   */
  cvec EvaluateH(rvec r);
  
  /**
   * @brief Integrate electric field against an RWG basis function.
   * @param f Pointer to RWG function.
   * @return Complex scalar integral value.
   */
  dcmplx IntegrateE(RWGFun* f);

  /**
   * @brief Integrate magnetic field against an RWG basis function.
   * @param f Pointer to RWG function.
   * @return Complex scalar integral value.
   */
  dcmplx IntegrateH(RWGFun* f);

  /**
   * @brief Integrate layered E and H fields against an RWG basis function.
   * @param f Pointer to RWG function.
   * @return Vector of complex values (E,H contributions).
   */
  std::vector<dcmplx> IntegrateLayeredEandH(RWGFun* f);
  
  /**
   * @brief Return the wavevector of the plane wave.
   * @return Complex wavevector.
   */
  cvec Wavevector();

  /**
   * @brief Return the polarization vector of the plane wave.
   * @return Polarization complex vector.
   */
  cvec Polarization();

  /**
   * @brief Write E and H field values at sample points to CSV files.
   * @param points List of 3D sample points.
   * @param eFilename Output filename for electric field.
   * @param hFilename Output filename for magnetic field.
   */
  void writeEHData(const std::vector<rvec>& points,
                   const std::string& eFilename,
                   const std::string& hFilename);

  /**
   * @brief Convert rvec to cvec (real vector to complex vector).
   * @param rv Real vector.
   * @return Complex vector.
   */
  cvec rvec2cvec(const rvec& rv) {
    return cvec(dcmplx(rv[0], 0.0), dcmplx(rv[1], 0.0), dcmplx(rv[2], 0.0));
  }
  
  /**
   * @brief Convert cvec to rvec, extracting real or imaginary parts.
   * @param cv Complex vector.
   * @param type "real" or "imag".
   * @return Real vector with corresponding components.
   */
  rvec cvec2rvec(const cvec& cv, std::string& type) {
    if (type == "real") {
      return rvec(cv[0].real(), cv[1].real(), cv[2].real());
    } else if (type == "imag") {
      return rvec(cv[0].imag(), cv[1].imag(), cv[2].imag());
    } else {
      std::cerr << "Error: Invalid type for cvec to rvec conversion. Use 'real' or 'imag'.\n";
      return rvec(0.0, 0.0, 0.0);  // Return a zero vector in case of error
    }
  }

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
   * @brief Identify field as a PlaneWave. 
   * @return True.
   */
  bool IsPlaneWave() { return true; }

  /** 
   * @brief Identify field as not a Dipole. 
   * @return False.
   */
  bool IsDipole() { return false; }

  /** 
   * @brief Identify field as not a Gaussian beam. 
   * @return False.
   */
  bool IsGaussian() { return false; }
  
  /**
   * @brief Provide layered-medium Green's function and its tabulation grids.
   * @param fGrnFields Pointer to a Green's-function backend (homogeneous or layered).
   * @param tabGrids Interpolation grids used by the layered backend.
   */
  virtual void setGrnFunLayered(GreenF *fGrnFields, std::vector<Grid> tabGrids);
};

#endif
