/**
 * @file GreenFHom3D.cpp
 * @brief Implementation of the 3D Green's Function class in a homogeneous background.
 */

#include "GreenFHom3D.h"
#include <iostream>
#include <tuple>
#include "quadrature.h"

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------
GreenFHom3D::GreenFHom3D(dcmplx wvl, dcmplx eps) : GreenF(wvl, eps) {}

GreenFHom3D::GreenFHom3D(dcmplx wvl, dcmplx eps, dcmplx mu)
    : GreenF(wvl, eps, mu) {}

GreenFHom3D::GreenFHom3D() : GreenF(dcmplx(0, 0), dcmplx(0, 0)) {}

// ---------------------------------------------------------------------------
// Proximity tests (heuristics to pick smoothed vs. raw kernels)
// ---------------------------------------------------------------------------
/**
 * @details Above \f$\sqrt{0.4}\f$ we are safely away from the \f$1/|r-r'|\f$ 
 * singular behavior, so the plain kernel is accurate enough without special care.
 */
bool GreenFHom3D::AboveThreshold(rvec r, rvec rp) {
  return norm(k_B) * dot(r - rp, r - rp) > 0.4;
}

/**
 * @details The distance between triangle centroids is compared to the average size (perimeter/2).
 * If the separation is larger than this size measure, simple quadrature is sufficient.
 */
bool GreenFHom3D::AboveThreshold(Triangle* T, Triangle* Tp) {
  double d = sqrt(dot(T->Center() - Tp->Center(), T->Center() - Tp->Center()));
  double psum = T->Perimeter() + Tp->Perimeter();
  return d > psum / 2.;
}

// ---------------------------------------------------------------------------
// Pointwise kernels
// ---------------------------------------------------------------------------
// Scalar Green's function G(r,r').
dcmplx GreenFHom3D::Evaluate(rvec r, rvec rp) {
  rvec dif(r - rp);
  double dis(sqrt(dot(dif, dif)));
  return exp(I * k_B * dis) / (4. * PI * dis);
}

// Dyadic Green's function \f$\overline{\overline G}(r,r')\f$.
cdyad GreenFHom3D::EvaluateDyadic(rvec r, rvec rp) {
  double dx(r(0) - rp(0));
  double dy(r(1) - rp(1));
  double dz(r(2) - rp(2));
  double rr(sqrt(dx * dx + dy * dy + dz * dz));
  dcmplx ff((I * k_B * rr - 1.) / (k_B * k_B * rr * rr));
  dcmplx nf((3. - 3. * I * k_B * rr - k_B * k_B * rr * rr) /
            (k_B * k_B * rr * rr * rr * rr));
  return cdyad(cvec(1. + ff + nf * dx * dx, nf * dx * dy, nf * dx * dz),
               cvec(nf * dy * dx, 1. + ff + nf * dy * dy, nf * dy * dz),
               cvec(nf * dz * dx, nf * dz * dy, 1. + ff + nf * dz * dz)) *
         exp(I * k_B * rr) / (4. * PI * rr);
}

/**
 * @details For small \f$|r-r'|\f$, we add/subtract a local expansion to tame the 
 * \f$1/|r-r'|\f$ singularity during area integrations. This returns 
 * \f$G(r,r') - 1/(4\pi|r-r'|) + k^2 |r-r'| / (8\pi)\f$ which is bounded at 
 * coalescence and integrates well over triangles.
 */
dcmplx GreenFHom3D::Smoothed(rvec r, rvec rp) {
  rvec dif(r - rp);
  double dis(sqrt(dot(dif, dif)));
  dcmplx out(I * k_B / (4. * PI));
  dcmplx id(1., 0.);
  if (dis != 0.)
    out = 1. / (4. * PI) *
          ((exp(I * k_B * dis) - id) / dis + k_B * k_B * dis / 2.);
  return out;
}

/**
 * @details Keeps only the \f$G(r,r') - 1/(4\pi|r-r'|)\f$ part. Used when 
 * triangles coincide (tristatus==3) to improve numerical stability without 
 * over-correcting the local behavior.
 */
dcmplx GreenFHom3D::halfSmoothed(rvec r, rvec rp) {
  rvec dif(r - rp);
  double dis(sqrt(dot(dif, dif)));
  dcmplx out(I * k_B / (4. * PI));
  dcmplx id(1., 0.);
  if (dis != 0.) out = 1. / (4. * PI) * ((exp(I * k_B * dis) - id) / dis);
  return out;
}

/**
 * @details Returns \f$\nabla^{\prime} G(r,r')\f$.
 */
cvec GreenFHom3D::Gradient(rvec r, rvec rp) {
  rvec dif(r - rp);
  double dis(sqrt(dot(dif, dif)));
  cvec out(-exp(I * k_B * dis) / (4. * PI * dis * dis) * (-1. / dis + I * k_B) *
           dif);
  return out;
}

/**
 * @details Useful to fill blocks that need G, gradG and dyadG consistently at the same point.
 */
std::tuple<dcmplx, cvec, cdyad> GreenFHom3D::EvaluateAll(rvec r, rvec rp) {
  rvec dif(r - rp);
  double dx(dif(0));
  double dy(dif(1));
  double dz(dif(2));
  double rr(sqrt(dx * dx + dy * dy + dz * dz));
  dcmplx t1(exp(I * k_B * rr) / (4. * PI * rr));
  dcmplx t2((I * k_B * rr - 1.) / (rr * rr));
  dcmplx t3((3. - 3. * I * k_B * rr - k_B * k_B * rr * rr) /
            (rr * rr * rr * rr));
  std::tuple<dcmplx, cvec, cdyad> ret;
  std::get<0>(ret) = t1;
  std::get<1>(ret) = -t1 * t2 * dif;
  std::get<2>(ret) =
      cdyad(cvec(t2 + t3 * dx * dx, t3 * dx * dy, t3 * dx * dz),
            cvec(t3 * dy * dx, t2 + t3 * dy * dy, t3 * dy * dz),
            cvec(t3 * dz * dx, t3 * dz * dy, t2 + t3 * dz * dz)) *
      t1;
  return ret;
}

// Smoothed gradient kernel.
cvec GreenFHom3D::GradientSmoothed(rvec r, rvec rp) {
  rvec dif(r - rp);
  double dis(sqrt(dot(dif, dif)));
  cvec out(0., 0., 0.);
  dcmplx id(1., 0.);
  if (dis != 0)
    out = -dif / (4. * PI * dis) *
          ((-exp(I * k_B * dis) + id) / (dis * dis) +
           I * k_B * exp(I * k_B * dis) / dis + k_B * k_B / 2.);
  return out;
}

// ---------------------------------------------------------------------------
// Triangle and post‑processing integrals
// ---------------------------------------------------------------------------
dcmplx GreenFHom3D::IntegrateD(RWGFun* /*f*/, RWGFun* /*fp*/) { exit(1); }
dcmplx GreenFHom3D::IntegrateK(RWGFun* /*f*/, RWGFun* /*fp*/) { exit(1); }

/**
 * @details Computes \f$\int_{T'} [ f_{p}(r') \; G(r,r') - \frac{\nabla' G(r,r')}{k^2}\,\mathrm{div} f_p ] \; dS' \f$,
 * choosing quadrature and (smoothed vs exact) kernels based on distance heuristics.
 * When near-singular, applies analytical singularity subtraction via \c SingSub.
 */
cvec GreenFHom3D::IntegrateDelta(rvec r, RWGFun* fp) {
  cvec out(0., 0., 0.);
  Triangle* Tp = fp->TrianglePtr();
  rvec rp[3];
  for (int i = 0; i < 3; i++) {
    rp[i] = Tp->Node(i);
  }
  double Ap = Tp->Area(), perimp = Tp->Perimeter();
  rvec Tc(Tp->Center());
  double dist = sqrt(dot(r - Tc, r - Tc));
  double kr = sqrt(norm(k_B)) * dist;

  int Ndunavant;
  double(*xdunavant)[3], *wdunavant;
  if (NeedsAccurate) {
    Ndunavant = (dist < perimp) ? dunavant::N17 : dunavant::N5;
    xdunavant = (dist < perimp) ? dunavant::x17 : dunavant::x5;
    wdunavant = (dist < perimp) ? dunavant::w17 : dunavant::w5;
  } else {
    Ndunavant = dunavant::N1;
    xdunavant = dunavant::x1;
    wdunavant = dunavant::w1;
  }

  bool NeedsSingSub = (dist < perimp) || (kr < 0.2 * PI);
  dcmplx (GreenFHom3D::*Dfunction)(rvec r, rvec rp) =
      NeedsSingSub ? &GreenFHom3D::Smoothed : &GreenFHom3D::Evaluate;
  cvec (GreenFHom3D::*Kfunction)(rvec r, rvec rp) =
      NeedsSingSub ? &GreenFHom3D::GradientSmoothed : &GreenFHom3D::Gradient;

  dcmplx Dtemp;
  cvec Ktemp;
  for (int i = 0; i < Ndunavant; ++i) {
    rvec sp = rp[0] * xdunavant[i][0] + rp[1] * xdunavant[i][1] +
              rp[2] * xdunavant[i][2];
    double dAp = Ap * wdunavant[i];
    Dtemp = dAp * (this->*Dfunction)(r, sp);
    Ktemp = dAp * (this->*Kfunction)(r, sp);
    rvec fpsp(fp->Evaluate(sp));
    out += (fpsp * Dtemp - 1. / (k_B * k_B) * fp->RWGPreFactor() * Ktemp);
  }

  if (NeedsSingSub) {
    SingSub ss;
    ss.SetTriangleRWG(fp, r);
    out += 1. / (4. * PI) * (-1. / (k_B * k_B) * fp->RWGPreFactor() *
                                 (ss.K3(-1) - k_B * k_B / 2. * ss.K3(1)) +
                             (ss.K2(-1) - k_B * k_B / 2. * ss.K2(1)));
  }

  return out;
}

/**
 * @details Computes \f$\int_{T'} f_p(r') \times \nabla' G(r,r')\, dS'\f$ with the same strategy
 * for quadrature / smoothing as in @ref IntegrateDelta, and uses the corresponding
 * K4 analytic terms for singularity subtraction.
 */
cvec GreenFHom3D::IntegrateKappa(rvec r, RWGFun* fp) {
  cvec out(0., 0., 0.);
  Triangle* Tp = fp->TrianglePtr();
  rvec rp[3];
  for (int i = 0; i < 3; i++) {
    rp[i] = Tp->Node(i);
  }
  double Ap = Tp->Area(), perimp = Tp->Perimeter();
  rvec Tc(Tp->Center());
  double dist = sqrt(dot(r - Tc, r - Tc));
  double kr = sqrt(norm(k_B)) * dist;

  int Ndunavant;
  double(*xdunavant)[3], *wdunavant;
  if (NeedsAccurate) {
    Ndunavant = (dist < perimp) ? dunavant::N17 : dunavant::N5;
    xdunavant = (dist < perimp) ? dunavant::x17 : dunavant::x5;
    wdunavant = (dist < perimp) ? dunavant::w17 : dunavant::w5;
  } else {
    Ndunavant = dunavant::N1;
    xdunavant = dunavant::x1;
    wdunavant = dunavant::w1;
  }

  bool NeedsSingSub = (dist < perimp) || (kr < 0.2 * PI);
  cvec (GreenFHom3D::*Kfunction)(rvec r, rvec rp) =
      NeedsSingSub ? &GreenFHom3D::GradientSmoothed : &GreenFHom3D::Gradient;

  dcmplx Dtemp;
  cvec Ktemp;
  for (int i = 0; i < Ndunavant; ++i) {
    rvec sp = rp[0] * xdunavant[i][0] + rp[1] * xdunavant[i][1] +
              rp[2] * xdunavant[i][2];
    double dAp = Ap * wdunavant[i];
    Ktemp = dAp * (this->*Kfunction)(r, sp);
    rvec fpsp(fp->Evaluate(sp));
    out += cross(Ktemp, (cvec)fpsp);
  }

  if (NeedsSingSub) {
    // Singularity subtraction for analytical integration on mesh triangles
    SingSub ss;
    ss.SetTriangleRWG(fp, r);
    out += 1. / (4. * PI) * (ss.K4(-1) - k_B * k_B / 2. * ss.K4(1));
  }

  return out;
}

// -----------------------------
// Pair utilities (matrix symmetry helpers)
// -----------------------------
// D-pair (symmetric for this homogeneous case).
std::pair<dcmplx, dcmplx> GreenFHom3D::PairD(RWGFun* f, RWGFun* fp) {
  // Don't care about symmetry, this is more exact!
  dcmplx elem(IntegrateD(f, fp));
  if (f == fp) return std::pair<dcmplx, dcmplx>(elem, elem);
  dcmplx elem2(IntegrateD(fp, f));
  return std::pair<dcmplx, dcmplx>(elem, elem2);
}

// K-pair (in general non-symmetric).
std::pair<dcmplx, dcmplx> GreenFHom3D::PairK(RWGFun* f, RWGFun* fp) {
  // Unfortunately, IntegrateK is not symmetric.
  dcmplx elem(IntegrateK(f, fp));
  if (f == fp) return std::pair<dcmplx, dcmplx>(elem, elem);
  dcmplx elem2(IntegrateK(fp, f));
  return std::pair<dcmplx, dcmplx>(elem, elem2);
}

// -----------------------------
// Assembly helpers
// -----------------------------
/**
 * @details This routine reuses geometry, quadrature points and kernel evaluations to fill many
 * matrix entries at once: \c fvec lives on triangle T, \c fpvec on triangle T'.
 * The logic mirrors the single-pair routines but amortizes the cost across pairs.
 */
void GreenFHom3D::SameTriDK(
    const std::vector<RWGFun*>& fvec, const std::vector<RWGFun*>& fpvec,
    std::vector<std::vector<std::pair<dcmplx, dcmplx>>>& DK) {
  Triangle *tPtr = fvec[0]->TrianglePtr(), *tpPtr = fpvec[0]->TrianglePtr();
  rvec r[3], rp[3];
  for (int i = 0; i < 3; i++) {
    r[i] = tPtr->Node(i);
    rp[i] = tpPtr->Node(i);
  }
  double A = tPtr->Area(), Ap = tpPtr->Area();
  int tristatus = tPtr->FindAdjacency(tpPtr);
  // Singularity subtraction for analytical integration on mesh triangles
  SingSub ss;

  int Ndunavant;
  double(*xdunavant)[3];
  double* wdunavant;
  bool NeedsDanalytic(false);
  bool NeedsKIntegral(tristatus != 3);
  bool NeedsKLineIntegral(false);
  bool NeedsSingSub;
  cvec (GreenFHom3D::*Kfunction)(rvec r, rvec rp);
  dcmplx (GreenFHom3D::*Dfunction)(rvec r, rvec rp);
  if (NeedsAccurate) {
    Ndunavant = tristatus > 0 ? dunavant::N17 : dunavant::N5;
    xdunavant = tristatus > 0 ? dunavant::x17 : dunavant::x5;
    wdunavant = tristatus > 0 ? dunavant::w17 : dunavant::w5;
    NeedsDanalytic = tristatus == 3;
    NeedsKLineIntegral = tristatus == 2;
    bool DistantAboveThreshold = tristatus == 0 and AboveThreshold(tPtr, tpPtr);
    NeedsSingSub = !(tristatus == 3 or DistantAboveThreshold);
    if (tristatus == 3) {
      Dfunction = &GreenFHom3D::halfSmoothed;
    } else if (!DistantAboveThreshold) {
      Dfunction = &GreenFHom3D::Smoothed;
    } else {
      Dfunction = &GreenFHom3D::Evaluate;
    }
    Kfunction = DistantAboveThreshold ? &GreenFHom3D::Gradient
                                      : &GreenFHom3D::GradientSmoothed;
  } else {
    Ndunavant = dunavant::N1;
    xdunavant = dunavant::x1;
    wdunavant = dunavant::w1;
    NeedsSingSub = !AboveThreshold(tPtr->Center(), tpPtr->Center());
    Kfunction =
        NeedsSingSub ? &GreenFHom3D::GradientSmoothed : &GreenFHom3D::Gradient;
    Dfunction = NeedsSingSub ? &GreenFHom3D::Smoothed : &GreenFHom3D::Evaluate;
  }

  std::vector<rvec> svec(Ndunavant), spvec(Ndunavant);
  for (int i = 0; i < Ndunavant; ++i) {
    svec[i] = r[0] * xdunavant[i][0] + r[1] * xdunavant[i][1] +
              r[2] * xdunavant[i][2];
    spvec[i] = rp[0] * xdunavant[i][0] + rp[1] * xdunavant[i][1] +
               rp[2] * xdunavant[i][2];
  }

  dcmplx Dtemp;
  cvec Ktemp;
  for (int i = 0; i < Ndunavant; ++i) {
    double dA = A * wdunavant[i];
    for (int j = 0; j < Ndunavant; ++j) {
      double dAp = Ap * wdunavant[j];
      Dtemp = dA * dAp * (this->*Dfunction)(svec[i], spvec[j]);
      if (NeedsKIntegral) {
        Ktemp = dA * dAp * (this->*Kfunction)(svec[i], spvec[j]);
      }
      for (int k = 0; k < (int)fvec.size(); ++k) {
        RWGFun* f = fvec[k];
        rvec fs(f->Evaluate(svec[i]));
        double divf(f->RWGPreFactor());
        for (int l = 0; l < (int)fpvec.size(); ++l) {
          RWGFun* fp = fpvec[l];
          rvec fpsp(fp->Evaluate(spvec[j]));
          double divfp(fp->RWGPreFactor());

          DK[k][l].first +=
              Dtemp * (dot(fs, fpsp) - 1. / (k_B * k_B) * divf * divfp);

          if (NeedsKIntegral) {
            DK[k][l].second += dot(fs, cross(Ktemp, (cvec)fpsp));
          }
        }
      }
    }

    if (NeedsSingSub) {
      ss.SetTriangleRWG(fpvec[0], svec[i]);
      for (int l = 0; l < (int)fpvec.size(); ++l) {
        RWGFun* fp = fpvec[l];
        double divfp(fp->RWGPreFactor());
        ss.ChangeTriangleRWG(fp);
        for (int k = 0; k < (int)fvec.size(); ++k) {
          RWGFun* f = fvec[k];
          double divf(f->RWGPreFactor());
          rvec fs(f->Evaluate(svec[i]));

          DK[k][l].first +=
              dA / (4. * PI) * (-1. / (k_B * k_B) * divf *
                                    (ss.K1(-1) - k_B * k_B / 2. * ss.K1(1)) +
                                dot(fs, ss.K2(-1) - k_B * k_B / 2. * ss.K2(1)));

          if (NeedsKLineIntegral) {
            DK[k][l].second +=
                dA / (4. * PI) *
                dot(fs, cross((rvec)(f->FreeVertex() - fp->FreeVertex()),
                              ss.K3n(-1)) *
                                divfp / (-2.) -
                            k_B * k_B / 2. * ss.K4(1));
          } else if (NeedsKIntegral) {
            DK[k][l].second +=
                dA / (4. * PI) * dot(fs, ss.K4(-1) - k_B * k_B / 2. * ss.K4(1));
          }
        }
      }
    }
  }

  if (NeedsKLineIntegral) {
    // Line integral for surface gradient of 1/R
    // Note that triangle order is reversed
    rvec n = cross((rvec)(rp[1] - rp[0]), (rvec)(rp[2] - rp[0]));
    for (int i = 0; i < 3; i++) {
      rvec &p1 = rp[i], &p2 = rp[(i + 1) % 3];
      double L = sqrt(dot(p2 - p1, p2 - p1));

      // Outward normal on edge;
      rvec m = cross((rvec)(p2 - p1), n);
      m = m / sqrt(dot(m, m));

      // Line integral
      for (int j = 0; j < gausslegendre::N17; j++) {
        rvec s = p1 * gausslegendre::x17[j][0] + p2 * gausslegendre::x17[j][1];
        double dL = L * gausslegendre::w17[j];

        ss.SetTriangleRWG(fvec[0], s);
        for (int k = 0; k < (int)fvec.size(); ++k) {
          RWGFun* f = fvec[k];
          ss.ChangeTriangleRWG(f);
          for (int l = 0; l < (int)fpvec.size(); ++l) {
            RWGFun* fp = fpvec[l];
            rvec pqm = cross((rvec)(f->FreeVertex() - fp->FreeVertex()), m);

            DK[k][l].second +=
                -dL / (4. * PI) * fp->RWGPreFactor() / 2. * dot(pqm, ss.K2(-1));
          }
        }
      }
    }
  }

  if (NeedsDanalytic) {
    for (int k = 0; k < (int)fvec.size(); ++k) {
      RWGFun* f = fvec[k];
      double divf(f->RWGPreFactor());
      for (int l = 0; l < (int)fpvec.size(); ++l) {
        RWGFun* fp = fpvec[l];
        double divfp(fp->RWGPreFactor());
        // Analytical 4D integral for (div f)(div f')/(4 pi r)
        DK[k][l].first +=
            -1. / (k_B * k_B) * divf * divfp / (4 * PI) * tPtr->getI1();
        // Analytical 4D integral for f dot f'/(4 pi r)
        // Note that f=fp and f!=fp need to be treated differently
        if (f == fp) {
          DK[k][l].first +=
              divf * divfp / (16. * PI) * tPtr->getI2(f->FreeVertexPtr());
        } else {
          DK[k][l].first +=
              divf * divfp / (16. * PI) *
              tPtr->getI2(f->FreeVertexPtr(), fp->FreeVertexPtr());
        }
      }
    }
  }
}

/**
 * @details DKG[k][l].first accumulates nablaG-weighted D-like terms; DKG[k][l].second holds
 * a 3-vector with row-wise contributions using the second-derivative dyad.
 */
void GreenFHom3D::SameTriDKG(
    const std::vector<RWGFun*>& fvec, const std::vector<RWGFun*>& fpvec,
    std::vector<std::vector<std::pair<dcmplx, dcmplx>>>& DK,
    std::vector<std::vector<std::pair<cvec, cvec>>>& DKG) {
  // Since DKG should be called only for different domains, we assume no need of
  // singularity subtraction!

  Triangle *tPtr = fvec[0]->TrianglePtr(), *tpPtr = fpvec[0]->TrianglePtr();
  rvec r[3], rp[3];
  for (int i = 0; i < 3; i++) {
    r[i] = tPtr->Node(i);
    rp[i] = tpPtr->Node(i);
  }
  double A = tPtr->Area(), Ap = tpPtr->Area();

  int Ndunavant;
  double(*xdunavant)[3];
  double* wdunavant;
  if (NeedsAccurate) {
    Ndunavant = dunavant::N5;
    xdunavant = dunavant::x5;
    wdunavant = dunavant::w5;
  } else {
    Ndunavant = dunavant::N1;
    xdunavant = dunavant::x1;
    wdunavant = dunavant::w1;
  }

  std::vector<rvec> svec(Ndunavant), spvec(Ndunavant);
  for (int i = 0; i < Ndunavant; ++i) {
    svec[i] = r[0] * xdunavant[i][0] + r[1] * xdunavant[i][1] +
              r[2] * xdunavant[i][2];
    spvec[i] = rp[0] * xdunavant[i][0] + rp[1] * xdunavant[i][1] +
               rp[2] * xdunavant[i][2];
  }

  std::tuple<dcmplx, cvec, cdyad> allG;
  for (int i = 0; i < Ndunavant; ++i) {
    double dA = A * wdunavant[i];
    for (int j = 0; j < Ndunavant; ++j) {
      double dAp = Ap * wdunavant[j];
      allG = EvaluateAll(svec[i], spvec[j]);
      std::get<0>(allG) *= dA * dAp;
      std::get<1>(allG) *= dA * dAp;
      std::get<2>(allG) *= dA * dAp;
      for (int k = 0; k < (int)fvec.size(); ++k) {
        RWGFun* f = fvec[k];
        rvec fs(f->Evaluate(svec[i]));
        double divf(f->RWGPreFactor());
        for (int l = 0; l < (int)fpvec.size(); ++l) {
          RWGFun* fp = fpvec[l];
          rvec fpsp(fp->Evaluate(spvec[j]));
          double divfp(fp->RWGPreFactor());

          DK[k][l].first += std::get<0>(allG) *
                            (dot(fs, fpsp) - 1. / (k_B * k_B) * divf * divfp);
          DK[k][l].second += dot(fs, cross(std::get<1>(allG), (cvec)fpsp));
          DKG[k][l].first += std::get<1>(allG) *
                             (dot(fs, fpsp) - 1. / (k_B * k_B) * divf * divfp);
          DKG[k][l].second[0] +=
              dot(fs, cross(std::get<2>(allG)[0], (cvec)fpsp));
          DKG[k][l].second[1] +=
              dot(fs, cross(std::get<2>(allG)[1], (cvec)fpsp));
          DKG[k][l].second[2] +=
              dot(fs, cross(std::get<2>(allG)[2], (cvec)fpsp));
        }
      }
    }
  }
}

//Compute delta and kappa vectors at an observation point for multiple RWGs on the same triangle.
void GreenFHom3D::SameTriDeltaKappa(
    rvec r, const std::vector<RWGFun*>& fpvec,
    std::vector<std::pair<cvec, cvec>>& DeltaKappa) {
  Triangle* Tp = fpvec[0]->TrianglePtr();
  rvec rp[3];
  for (int i = 0; i < 3; i++) {
    rp[i] = Tp->Node(i);
  }
  double Ap = Tp->Area(), perimp = Tp->Perimeter();
  rvec Tc(Tp->Center());
  double dist = sqrt(dot(r - Tc, r - Tc));
  double kr = sqrt(norm(k_B)) * dist;
  SingSub ss;

  int Ndunavant;
  double(*xdunavant)[3], *wdunavant;
  if (NeedsAccurate) {
    Ndunavant = (dist < perimp) ? dunavant::N17 : dunavant::N5;
    xdunavant = (dist < perimp) ? dunavant::x17 : dunavant::x5;
    wdunavant = (dist < perimp) ? dunavant::w17 : dunavant::w5;
  } else {
    Ndunavant = dunavant::N1;
    xdunavant = dunavant::x1;
    wdunavant = dunavant::w1;
  }

  bool NeedsSingSub = (dist < perimp) || (kr < 0.2 * PI);
  dcmplx (GreenFHom3D::*Dfunction)(rvec r, rvec rp) =
      NeedsSingSub ? &GreenFHom3D::Smoothed : &GreenFHom3D::Evaluate;
  cvec (GreenFHom3D::*Kfunction)(rvec r, rvec rp) =
      NeedsSingSub ? &GreenFHom3D::GradientSmoothed : &GreenFHom3D::Gradient;

  dcmplx Dtemp;
  cvec Ktemp;
  for (int i = 0; i < Ndunavant; ++i) {
    rvec sp = rp[0] * xdunavant[i][0] + rp[1] * xdunavant[i][1] +
              rp[2] * xdunavant[i][2];
    double dAp = Ap * wdunavant[i];
    Dtemp = dAp * (this->*Dfunction)(r, sp);
    Ktemp = dAp * (this->*Kfunction)(r, sp);
    for (int j = 0; j < (int)fpvec.size(); ++j) {
      RWGFun* fp = fpvec[j];
      rvec fpsp(fp->Evaluate(sp));
      DeltaKappa[j].first +=
          (fpsp * Dtemp - 1. / (k_B * k_B) * fp->RWGPreFactor() * Ktemp);
      DeltaKappa[j].second += cross(Ktemp, (cvec)fpsp);
    }
  }

  if (NeedsSingSub) {
    ss.SetTriangleRWG(fpvec[0], r);
    for (int j = 0; j < (int)fpvec.size(); ++j) {
      RWGFun* fp = fpvec[j];
      ss.ChangeTriangleRWG(fp);
      DeltaKappa[j].first +=
          1. / (4. * PI) * (-1. / (k_B * k_B) * fp->RWGPreFactor() *
                                (ss.K3(-1) - k_B * k_B / 2. * ss.K3(1)) +
                            (ss.K2(-1) - k_B * k_B / 2. * ss.K2(1)));
      DeltaKappa[j].second +=
          1. / (4. * PI) * (ss.K4(-1) - k_B * k_B / 2. * ss.K4(1));
    }
  }
}
