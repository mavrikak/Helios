/**
 * @file WaveguideMode.cpp
 * @brief Implementation for the WaveguideMode class.
 */

#include "WaveguideMode.h"
#include <stdexcept>
#include <algorithm>
#include <iostream>

WaveguideMode::WaveguideMode(dcmplx wvl, std::vector<dcmplx> k_L, std::vector<dcmplx> epsL,
                             std::vector<dcmplx> muL, std::vector<double> zValsInterfaces,
                             std::vector<double> thickness, std::string polarization,
                             double propagation_angle_phi, int mode_nr)
    : wvl(wvl), pol(polarization), phi(propagation_angle_phi) {

    k0 = 2.0 * PI / wvl;
    layeredUtils = new LayeredMediaUtils(k_L, epsL, muL, zValsInterfaces, thickness);

    // Establish the valid search bounds for bound modes
    double kmin = 0, kmax = 0;
    for (const auto& k : k_L) {
        kmax = std::max(kmax, std::real(k));
    }
    kmin = std::max(std::real(k_L.front()), std::real(k_L.back()));

    double smin = kmin + 1e-9;
    double smax = kmax - 1e-9;

    if (smax <= smin) {
        throw std::runtime_error("WaveguideMode: No guided modes possible (core index not higher than cladding).");
    }

    // Scan and find the specified mode using bisection
    int mode_cnt = 0;
    int steps = 1000;
    double dk = (smax - smin) / steps;
    double xup = smax;
    dcmplx det_prev = detfunc(xup);
    double found_krho = -1.0;

    for (int i = 1; i <= steps; ++i) {
        double xlow = smax - i * dk;
        dcmplx det_curr = detfunc(xlow);

        if (std::real(det_prev) * std::real(det_curr) <= 0.0) {
            double r = bisection(xlow, xup, true);
            if (std::abs(std::imag(detfunc(r))) < std::abs(std::real(detfunc(r))) || std::abs(std::imag(detfunc(r))) < 1e-5) {
                if (mode_cnt == mode_nr) { found_krho = r; break; }
                mode_cnt++;
            }
        }
        if (found_krho < 0.0 && std::imag(det_prev) * std::imag(det_curr) <= 0.0) {
            double r = bisection(xlow, xup, false);
            if (std::abs(std::real(detfunc(r))) < std::abs(std::imag(detfunc(r))) || std::abs(std::real(detfunc(r))) < 1e-5) {
                if (mode_cnt == mode_nr) { found_krho = r; break; }
                mode_cnt++;
            }
        }
        if (found_krho > 0.0) break;
        
        xup = xlow;
        det_prev = det_curr;
    }

    if (found_krho < 0.0) {
        throw std::runtime_error("WaveguideMode: Requested mode index not found.");
    }
    krho = found_krho;

    // Compute global field amplitudes (A, B) for each layer using the transfer matrix method
    int n_layers = k_L.size();
    globalA.resize(n_layers);
    globalB.resize(n_layers);

    // Bound mode definition: exponential decay in the top cladding
    std::vector<dcmplx> V_curr = {1.0, 0.0}; 
    int top_layer = n_layers - 1;

    dcmplx kz_top = this->complexSqrt(k_L[top_layer] * k_L[top_layer] - krho * krho);
    globalA[top_layer] = V_curr[0] * std::exp(-I * kz_top * zValsInterfaces.back());
    globalB[top_layer] = V_curr[1] * std::exp(I * kz_top * zValsInterfaces.back());

    // Cascade down through interior layers
    for (int i = top_layer - 1; i >= 1; --i) {
        auto M = this->layeredUtils->InterfaceTransferMatrix(i, i + 1, krho, pol);
        std::vector<dcmplx> V_next(2);
        V_next[0] = M[0][0] * V_curr[0] + M[0][1] * V_curr[1];
        V_next[1] = M[1][0] * V_curr[0] + M[1][1] * V_curr[1];

        dcmplx kz_i = this->complexSqrt(k_L[i] * k_L[i] - krho * krho);
        globalA[i] = V_next[0] * std::exp(-I * kz_i * zValsInterfaces[i]);
        globalB[i] = V_next[1] * std::exp(I * kz_i * zValsInterfaces[i]);

        auto P = this->layeredUtils->propagationMatrix(i, krho);
        V_curr[0] = P[0][0] * V_next[0] + P[0][1] * V_next[1];
        V_curr[1] = P[1][0] * V_next[0] + P[1][1] * V_next[1];
    }

    // Connect to the bottom cladding
    auto M0 = this->layeredUtils->InterfaceTransferMatrix(0, 1, krho, pol);
    std::vector<dcmplx> V_0(2);
    V_0[0] = M0[0][0] * V_curr[0] + M0[0][1] * V_curr[1];
    V_0[1] = M0[1][0] * V_curr[0] + M0[1][1] * V_curr[1];

    dcmplx kz_0 = this->complexSqrt(k_L[0] * k_L[0] - krho * krho);
    globalA[0] = V_0[0] * std::exp(-I * kz_0 * zValsInterfaces[0]);
    globalB[0] = V_0[1] * std::exp(I * kz_0 * zValsInterfaces[0]);

    // Normalize the mode power
    // Calculate the evanescent decay lengths to bound the integration window
    double im_bot = std::abs(std::imag(kz_0));
    double im_top = std::abs(std::imag(kz_top));
    double decay_bot = (im_bot > 1e-6) ? 1.0 / im_bot : 1e6;
    double decay_top = (im_top > 1e-6) ? 1.0 / im_top : 1e6;

    double z_start = zValsInterfaces.front() - 5.0 * decay_bot;
    double z_end = zValsInterfaces.back() + 5.0 * decay_top;
    
    int n_pts = 10000;
    double dz = (z_end - z_start) / n_pts;
    double P_total = 0.0;

    // Integrate the Poynting vector along the direction of propagation using the trapezoidal rule
    for (int i = 0; i <= n_pts; ++i) {
        double z = z_start + i * dz;
        rvec r(0.0, 0.0, z); 
        auto fields = EvaluateLayeredEandH(r);
        cvec E = fields[0];
        cvec H = fields[1];

        dcmplx Sx = 0.5 * (E(1) * std::conj(H(2)) - E(2) * std::conj(H(1)));
        dcmplx Sy = 0.5 * (E(2) * std::conj(H(0)) - E(0) * std::conj(H(2)));
        double S_rho = std::real(Sx) * std::cos(phi) + std::real(Sy) * std::sin(phi);

        double weight = (i == 0 || i == n_pts) ? 0.5 : 1.0;
        P_total += S_rho * weight * dz;
    }

    // Scale all layer amplitudes by the normalization factor
    double norm_factor = 1.0 / std::sqrt(std::abs(P_total));
    for (int i = 0; i < n_layers; ++i) {
        globalA[i] *= norm_factor;
        globalB[i] *= norm_factor;
    }    
}

dcmplx WaveguideMode::detfunc(double krho_test) {
    auto mTot = this->layeredUtils->totalTransferMatrix(krho_test, pol);
    return mTot[0][0];
}

double WaveguideMode::bisection(double a, double b, bool use_real) {
    double fa = use_real ? std::real(detfunc(a)) : std::imag(detfunc(a));
    for (int i = 0; i < 100; ++i) {
        double c = (a + b) / 2.0;
        double fc = use_real ? std::real(detfunc(c)) : std::imag(detfunc(c));
        if (std::abs(fc) < 1e-12 || (b - a) < 1e-12) return c;
        if (fa * fc <= 0.0) { b = c; }
        else { a = c; fa = fc; }
    }
    return (a + b) / 2.0;
}

cvec WaveguideMode::EvaluateE(rvec r) {
    return EvaluateLayeredEandH(r)[0];
}

cvec WaveguideMode::EvaluateH(rvec r) {
    return EvaluateLayeredEandH(r)[1];
}

std::vector<cvec> WaveguideMode::EvaluateLayeredEandH(rvec r) {
    int i = this->layeredUtils->GetLayerIndex(r);
    dcmplx kz = this->complexSqrt(this->layeredUtils->getkL(i) * this->layeredUtils->getkL(i) - krho * krho);
    dcmplx A = globalA[i];
    dcmplx B = globalB[i];

    dcmplx phaseA = std::exp(I * kz * r[2]);
    dcmplx phaseB = std::exp(-I * kz * r[2]);
    dcmplx xyPhase = std::exp(I * krho * (r[0] * std::cos(phi) + r[1] * std::sin(phi)));

    dcmplx E_x = 0, E_y = 0, E_z = 0;
    dcmplx H_x = 0, H_y = 0, H_z = 0;

    dcmplx eps_i = this->layeredUtils->getEps(i);
    dcmplx mu_i = this->layeredUtils->getMu(i);

    if (pol == "TE") {
        dcmplx Ey_A = A * phaseA;
        dcmplx Ey_B = B * phaseB;

        E_x = -std::sin(phi) * (Ey_A + Ey_B);
        E_y =  std::cos(phi) * (Ey_A + Ey_B);
        E_z = 0.0;

        dcmplx factor = 1.0 / (k0 * Z0 * mu_i);
        H_x = factor * (-kz * std::cos(phi) * Ey_A + kz * std::cos(phi) * Ey_B);
        H_y = factor * (-kz * std::sin(phi) * Ey_A + kz * std::sin(phi) * Ey_B);
        H_z = factor * (krho * (Ey_A + Ey_B));
    } else if (pol == "TM") {
        dcmplx Hy_A = A * phaseA;
        dcmplx Hy_B = B * phaseB;

        H_x = -std::sin(phi) * (Hy_A + Hy_B);
        H_y =  std::cos(phi) * (Hy_A + Hy_B);
        H_z = 0.0;

        dcmplx factor = Z0 / (k0 * eps_i);
        E_x = factor * (kz * std::cos(phi) * Hy_A - kz * std::cos(phi) * Hy_B);
        E_y = factor * (kz * std::sin(phi) * Hy_A - kz * std::sin(phi) * Hy_B);
        E_z = factor * (-krho * (Hy_A + Hy_B));
    }

    E_x *= xyPhase; E_y *= xyPhase; E_z *= xyPhase;
    H_x *= xyPhase; H_y *= xyPhase; H_z *= xyPhase;

    cvec E(E_x, E_y, E_z);
    cvec H(H_x, H_y, H_z);
    return {E, H};
}

dcmplx WaveguideMode::IntegrateE(RWGFun* f) {
    cvec (WaveguideMode::*EvalFunc)(rvec) = &WaveguideMode::EvaluateE;
    return gauss.Integrate(this, EvalFunc, f);
}

dcmplx WaveguideMode::IntegrateH(RWGFun* f) {
    cvec (WaveguideMode::*EvalFunc)(rvec) = &WaveguideMode::EvaluateH;
    return gauss.Integrate(this, EvalFunc, f);
}

std::vector<dcmplx> WaveguideMode::IntegrateLayeredEandH(RWGFun* f) {
    std::vector<cvec> (WaveguideMode::*EvalFunc)(rvec) = &WaveguideMode::EvaluateLayeredEandH;
    return gauss.Integrate(this, EvalFunc, f, NeedsAccurate);
}

void WaveguideMode::setGrnFunLayered(GreenF* fGrnFields, std::vector<Grid> tabGrids) {
    // Left empty: Guided modes directly yield the field and bypass the need
    // to trigger secondary tabular scattering evaluations.
}