/**
 * @file PlaneWave.cpp
 * @brief Implementation for the PlaneWave class.
 */

#include "PlaneWave.h"
#include <iostream>

/**
 * @details Normalizes the input propagation direction, sets the wavenumber
 * \f$\mathbf{k} = \frac{2\pi}{\lambda}\sqrt{\varepsilon_r}\,\hat{\mathbf{d}}\f$,
 * the medium impedance \f$ z = Z_0/\sqrt{\varepsilon_r}\f$, and stores polarization.
 */
PlaneWave::PlaneWave(dcmplx wvl, double eps, rvec dir, cvec po) {
  layeredUtils = nullptr; // No layered media utilities for this constructor
  dir /= sqrt(dot(dir, dir));
  k = (2 * PI * sqrt(eps) / wvl) * dir;
  z = Z0 / sqrt(eps);
  p = po;
}

/**
 * @details Allocates a LayeredMediaUtils helper and initializes propagation parameters
 * depending on whether the wave impinges from +z (top half-space) or −z (bottom).
 * Also computes the transverse wavenumber magnitude \f$k_\rho\f$.
 */
PlaneWave::PlaneWave(std::vector<dcmplx> k_L, std::vector<dcmplx> epsL,
                     std::vector<dcmplx> muL, std::vector<double> zValsInterfaces,
                     std::vector<double> thickness, rvec propDir, cvec po) {
  layeredUtils = new LayeredMediaUtils(k_L, epsL, muL, zValsInterfaces, thickness);
  propDir /= sqrt(dot(propDir, propDir));
  dir = propDir;
  z = 1.0;
  if (dir[2] > 0.0) {
    // Propagation in the positive z-direction
    k  = k_L[0] * dir;
    zL = Z0 * sqrt(muL[0] / epsL[0]);
    incLayerZ = zValsInterfaces[0];
    incLayerIndex = 0;
  } else {
    // Propagation in the negative z-direction
    k  = k_L.back() * dir;
    zL = Z0 * sqrt(muL.back() / epsL.back());
    incLayerZ = zValsInterfaces.back();
    incLayerIndex = zValsInterfaces.size();
  }
  p = po;
  krho = sqrt(k[0] * k[0] + k[1] * k[1]);
#ifdef DEBUG
  std::string filename = "pos_inc.txt";
  std::vector<rvec> points = this->layeredUtils->readPositionsFromFile(filename);
  writeEHData(points, "E_data.csv", "H_data.csv");
  for (int i = 0; i < points.size(); ++i) {
    std::vector<cvec> EandH = EvaluateLayeredEandH(points[i]);
    std::cout << "E" << i << " at point " << points[i] << ": " << EandH[0] << std::endl;
    std::cout << "H" << i << " at point " << points[i] << ": " << EandH[1] << std::endl;
  }
#endif
}

/**
 * @details Uses \f$\mathbf{E}(\mathbf{r}) = \mathbf{p}\,e^{j\mathbf{k}\cdot\mathbf{r}}\f$.
 */
cvec PlaneWave::EvaluateE(rvec r) { return p * exp(I * dot(k, r)); }

/**
 * @details Decomposes the incident polarization into local TE/TM unit vectors,
 * computes incident amplitudes, queries LayeredMediaUtils::Secondary
 * to obtain up/down "secondary" wave amplitudes within the observation layer,
 * and assembles the fields with impedance scaling in that layer. If the observation 
 * layer equals the incidence layer, the direct incident contribution is added on top 
 * of the secondary fields.
 */
std::vector<cvec> PlaneWave::EvaluateLayeredEandH(rvec r) {
  bool normalInc = std::abs(dir[2]) == 1.0;
  cvec te(dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0));
  cvec tm(dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0));
  cvec Einc(dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0));
  cvec Hinc(dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0));

  // Build TE/TM basis
  if (normalInc) {
    te = p / sqrt(dot(p, p));
    tm = cross(te, rvec2cvec(dir));
  } else {
    rvec temp = cross(dir, rvec(0.0, 0.0, 1.0));
    temp = temp / sqrt(dot(temp, temp));
    te = rvec2cvec(temp);
    tm = cross(te, rvec2cvec(dir));
  }

  // Incident amplitudes at the incidence interface
  std::map<std::string, dcmplx> amplitude;
  amplitude["TE"] = dot(te, p) * exp(I * k[2] * incLayerZ);
  amplitude["TM"] = dot(tm, p) * exp(I * k[2] * incLayerZ) / zL;
  
  // Observation layer and secondary coefficients
  int obsLayerIndex = this->layeredUtils->GetLayerIndex(r);
  std::map<std::string, blitz::Array<dcmplx, 2>> secondary;
  secondary.emplace("TE", this->layeredUtils->Secondary(krho, obsLayerIndex, incLayerIndex, "TE"));
  secondary.emplace("TM", this->layeredUtils->Secondary(krho, obsLayerIndex, incLayerIndex, "TM"));
  
  // Longitudinal component in the observation layer
  dcmplx kZ = this->layeredUtils->getkL(obsLayerIndex);
  kZ = this->complexSqrt(kZ * kZ - krho * krho);

  // Assemble E and H from secondary waves
  std::array<std::string,2> polarizations = { "TE", "TM" };
  cvec E(dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0));
  cvec H(dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0));

  // Sum contributions for TE and TM
  for (const auto& name : polarizations) {
    auto sec = secondary[name];
    for (int i = 0; i < sec.extent(0); ++i) {
      double dir1  = (i == 0 && obsLayerIndex != 0) ? 1.0 : -1.0;
      double zeta  = (i == 0 && obsLayerIndex != 0)
                   ? (r[2] - this->layeredUtils->getZValInterface(obsLayerIndex - 1))
                   : (this->layeredUtils->getZValInterface(obsLayerIndex) - r[2]);
      dcmplx xpar = k[0] * r[0] + k[1] * r[1] + kZ * zeta;
      cvec dirNew(k[0] / this->layeredUtils->getkL(obsLayerIndex),
                  k[1] / this->layeredUtils->getkL(obsLayerIndex),
                  dir1 * kZ / this->layeredUtils->getkL(obsLayerIndex));
      if (name == "TE") {
        E = te * amplitude[name] * sec(i, 0) * exp(I * xpar);
        H = cross(dirNew, E) / ( Z0 * sqrt( this->layeredUtils->getMu(obsLayerIndex) / 
                                            this->layeredUtils->getEps(obsLayerIndex) ) );
      } else if (name == "TM") {
        H = amplitude[name] * sec(i, 0) * exp(I * xpar) * cross(dir, p);
        E = -cross(dirNew, H) * ( Z0 * sqrt( this->layeredUtils->getMu(obsLayerIndex) / 
                                             this->layeredUtils->getEps(obsLayerIndex) ) );
      }
      Einc += E; Hinc += H;
    }
  }
  
  // Add direct incident field when observation == incidence layer
  if (obsLayerIndex == incLayerIndex) {
    Einc += this->EvaluateE(r);
    Hinc += this->EvaluateH(r) / 
            ( Z0 * sqrt( this->layeredUtils->getMu(obsLayerIndex) / 
                         this->layeredUtils->getEps(obsLayerIndex) ) );
  }
  return {Einc, Hinc};
}

/**
 * @details Uses \f$\mathbf{H}(\mathbf{r}) = \frac{1}{z}\,\hat{\mathbf{k}}\times\mathbf{p}\,e^{j\mathbf{k}\cdot\mathbf{r}}\f$.
 */
cvec PlaneWave::EvaluateH(rvec r) {
  cvec prop(k);
  prop /= sqrt(dot(prop, prop));
  cvec hp(cross(prop, p));
  return hp * exp(I * dot(k, r)) / z;
}

/**
 * @details Uses Gaussian quadrature through GaussQuad<PlaneWave>::Integrate.
 */
dcmplx PlaneWave::IntegrateE(RWGFun* f) {
  dcmplx out;
  cvec (PlaneWave::*EvaluationFunction)(rvec r);
  EvaluationFunction = &PlaneWave::EvaluateE;
  out = gauss.Integrate(this, EvaluationFunction, f);
  return out;
}


/**
 * @details Uses Gaussian quadrature through GaussQuad<PlaneWave>::Integrate.
 */
dcmplx PlaneWave::IntegrateH(RWGFun* f) {
  dcmplx out;
  cvec (PlaneWave::*EvaluationFunction)(rvec r);
  EvaluationFunction = &PlaneWave::EvaluateH;
  out = gauss.Integrate(this, EvaluationFunction, f);
  return out;
}

/**
 * @details Selects #EvaluateLayeredEandH for the integrand and passes the 
 * global #NeedsAccurate flag to the quadrature object.
 */
std::vector<dcmplx> PlaneWave::IntegrateLayeredEandH(RWGFun* f) {
  std::vector<dcmplx> out;
  std::vector<cvec> (PlaneWave::*EvaluationFunction)(rvec r);
  EvaluationFunction = &PlaneWave::EvaluateLayeredEandH;
  out = gauss.Integrate(this, EvaluationFunction, f, NeedsAccurate);
  return out;
}

// Get the wavevector \f$\mathbf{k}\f$.
cvec PlaneWave::Wavevector() { return k; }

// @brief Get the polarization vector \f$\mathbf{p}\f$.
cvec PlaneWave::Polarization() { return p; }

/**
 * @throws std::runtime_error always.
 */
void PlaneWave::setGrnFunLayered(GreenF* /*fGrnFields*/, std::vector<Grid> /*tabGrids*/) {
  throw std::runtime_error("PlaneWave does not support layered Green's functions.");
}

/**
 * @details Produces two files with headers:
 *  - E CSV: x,y,z,Ex_re,Ex_im,Ey_re,Ey_im,Ez_re,Ez_im
 *  - H CSV: x,y,z,Hx_re,Hx_im,Hy_re,Hy_im,Hz_re,Hz_im
 * @throws std::runtime_error if files cannot be opened for writing.
 */
void PlaneWave::writeEHData(
    const std::vector<rvec>& points,
    const std::string& eFilename,
    const std::string& hFilename)
{
    // 1) open output streams
    std::ofstream eFile(eFilename);
    std::ofstream hFile(hFilename);
    if (!eFile.is_open() || !hFile.is_open()) {
        throw std::runtime_error("Failed to open E/H output files");
    }

    // 2) write CSV headers
    eFile << "x,y,z,Ex_re,Ex_im,Ey_re,Ey_im,Ez_re,Ez_im\n";
    hFile << "x,y,z,Hx_re,Hx_im,Hy_re,Hy_im,Hz_re,Hz_im\n";

    // 3) loop over points
    for (size_t i = 0; i < points.size(); ++i) {
        const auto& pt = points[i];
        auto fields  = EvaluateLayeredEandH(pt);
        const cvec& E = fields[0];
        const cvec& H = fields[1];

        // write E line
        eFile
          << pt[0] << ',' << pt[1] << ',' << pt[2] << ','
          << std::real(E[0]) << ',' << std::imag(E[0]) << ','
          << std::real(E[1]) << ',' << std::imag(E[1]) << ','
          << std::real(E[2]) << ',' << std::imag(E[2]) << '\n';

        // write H line
        hFile
          << pt[0] << ',' << pt[1] << ',' << pt[2] << ','
          << std::real(H[0]) << ',' << std::imag(H[0]) << ','
          << std::real(H[1]) << ',' << std::imag(H[1]) << ','
          << std::real(H[2]) << ',' << std::imag(H[2]) << '\n';
    }

    // 4) close files
    eFile.close();
    hFile.close();
}
