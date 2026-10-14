/**
 * @file WaveguideMode.h
 * @brief Defines the WaveguideMode class: representation of a guided mode excitation.
 * Supports layered media. Solves for the propagation constant using the transfer matrix
 * method and evaluates the corresponding modal fields.
 */

#ifndef WAVEGUIDEMODE_H
#define WAVEGUIDEMODE_H

#include "GaussQuad.h"
#include "IncidentField.h"
#include "globals.h"
#include "LayeredMediaUtils.h"
#include <vector>
#include <string>

/// \ingroup incidentField
/**
 * @class WaveguideMode
 * @brief Models a guided waveguide mode as an incident field in a layered medium.
 */
class WaveguideMode : public IncidentField {
private:
    dcmplx wvl;                       ///< Incident wavelength
    dcmplx k0;                        ///< Vacuum wavenumber
    double krho;                      ///< Transverse wavenumber of the mode
    std::string pol;                  ///< Polarization ("TE" or "TM")
    double phi;                       ///< Propagation angle in the xy-plane
    std::vector<dcmplx> globalA;      ///< Forward-propagating amplitude in each layer
    std::vector<dcmplx> globalB;      ///< Backward-propagating amplitude in each layer
    LayeredMediaUtils* layeredUtils;  ///< Utilities for layered case
    GaussQuad<WaveguideMode> gauss;   ///< Gaussian quadrature integrator

    /**
     * @brief Evaluates the M_00 component of the total transfer matrix.
     * @param krho_test Test transverse wavenumber.
     * @return Complex M_00 value.
     */
    dcmplx detfunc(double krho_test);

    /**
     * @brief Bisection root-finding method for the dispersion relation.
     * @param a Lower bound of krho.
     * @param b Upper bound of krho.
     * @param use_real Whether to track the zero-crossing of the real or imaginary part.
     * @return The isolated root krho.
     */
    double bisection(double a, double b, bool use_real);

public:
    /**
     * @brief Construct a guided mode excitation.
     * @param wvl Vacuum wavelength.
     * @param k_L Wave numbers per layer.
     * @param epsL Permittivity values per layer.
     * @param muL Permeability values per layer.
     * @param zValsInterfaces z-coordinates of interfaces.
     * @param thickness Thicknesses of layers.
     * @param polarization "TE" or "TM".
     * @param propagation_angle_phi Angle of propagation relative to x-axis.
     * @param mode_nr The index of the guided mode to excite (0 = fundamental).
     */
    WaveguideMode(dcmplx wvl, std::vector<dcmplx> k_L, std::vector<dcmplx> epsL,
                  std::vector<dcmplx> muL, std::vector<double> zValsInterfaces,
                  std::vector<double> thickness, std::string polarization,
                  double propagation_angle_phi, int mode_nr);

    ~WaveguideMode() override {
        if (layeredUtils) { delete layeredUtils; }
    }

    cvec EvaluateE(rvec r) override;
    cvec EvaluateH(rvec r) override;
    std::vector<cvec> EvaluateLayeredEandH(rvec r) override;

    dcmplx IntegrateE(RWGFun* f) override;
    dcmplx IntegrateH(RWGFun* f) override;
    std::vector<dcmplx> IntegrateLayeredEandH(RWGFun* f) override;

    bool IsPlaneWave() override { return false; }
    bool IsDipole() override { return false; }
    bool IsGaussian() override { return false; }
    bool IsWaveguideMode() override { return true; }
    void setGrnFunLayered(GreenF* fGrnFields, std::vector<Grid> tabGrids) override;

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
};

#endif