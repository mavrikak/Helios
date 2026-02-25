/**
 * @file SIEFormPMCHW.cpp
 * @brief Implementation for the SIEFormPMCHW class.
 */

#include "SIEFormPMCHW.h"
#include <mutex>
#include <thread>

// ---------------- Constructors ----------------
// Construct PMCHWT for a homogeneous medium (no layered back-end).
SIEFormPMCHW::SIEFormPMCHW(std::vector<RWGFun> *inRWGFuns, GreenF *inGrnFun,
                           std::vector<IncidentField *> *inIncFields,
                           dcmplx inVacWavelength, dcmplx inEpsilon,
                           dcmplx inMu)
    : SIEForm(inRWGFuns, inGrnFun, inIncFields, inVacWavelength),
      epsilon(inEpsilon),
      mu(inMu) { grnFunLayered = nullptr; }


// Construct PMCHWT for a layered medium (requires layered Green's).
SIEFormPMCHW::SIEFormPMCHW(std::vector<RWGFun> *inRWGFuns,
                           GreenF *inGrnFun,
                           std::vector<IncidentField *> *inIncFields,
                           dcmplx inVacWavelength,
                           GreenF *inGrnFunLayered,
                           const std::vector<dcmplx>& inEpsilon,
                           const std::vector<dcmplx>& inMu,
                           const std::vector<double>& inZVals,
                           int inMidLayerIndex)
    : SIEForm(inRWGFuns, inGrnFun, inIncFields, inVacWavelength),
      epsilonLayered(inEpsilon),
      muLayered(inMu),
      zValsInterfaces(inZVals),
      grnFunLayered(inGrnFunLayered),
      midLayerIndex(inMidLayerIndex) {}

// ---------------- Medium helpers -----------------
// Homogeneous wave impedance Zr = sqrt(mu/epsilon).
dcmplx SIEFormPMCHW::Zr(dcmplx e, dcmplx m) { return sqrt(m / e); }

// Layer-wise impedances Zr_l = sqrt(μ_l/ε_l) with safe length.
std::vector<dcmplx> SIEFormPMCHW::ZrLayered(const std::vector<dcmplx>& e, 
                                            const std::vector<dcmplx>& m) {
  // Ensure both vectors have the same size
  size_t size = std::min(e.size(), m.size());
  std::vector<dcmplx> Zr_values(size);
  for (size_t i = 0; i < size; ++i) {
    Zr_values[i] = csqrt(m[i] / e[i]);  // Zr = sqrt(mu / epsilon)
  }
  return Zr_values;
}

// ---------------- Matrix assembly ----------------
/**
 * @details Partitions RWGs by owner triangle, spawns (threads-1) 
 * workers plus the current thread to call FillLSEMatrixParallel. 
 * Uses per-row mutexes to update matrix blocks safely.
 */
int SIEFormPMCHW::FillLSEMatrix(blitz::Array<dcmplx, 2> *lseMatrix) {
  using std::vector;
  using std::pair;
  vector<RWGFun>::iterator rwgM, rwgN;

  // Partition the RWG functions to groups over same triangle
  vector<vector<RWGFun *>> rwgPerTri;
  for (rwgM = rwgFuns->begin(); rwgM != rwgFuns->end();) {
    rwgPerTri.push_back(vector<RWGFun *>());
    for (rwgN = rwgM;
         rwgN != rwgFuns->end() and rwgN->TrianglePtr() == rwgM->TrianglePtr();
         rwgN++) {
      rwgPerTri.back().push_back(&(*rwgN));
    }
    rwgM = rwgN;
  }

#ifdef DEBUG
  for (int i = 0; i < rwgTriLen; ++i) {
    std::cout << "Group " << i << " : size " << rwgPerTri[i].size() << "\n";
  }
#endif

  std::cout << "Filling matrix with " << rwgFuns->size() << " elements."
            << std::endl;

  int mxOffset(lseMatrix->shape()(0) / 2);
  std::mutex *mut = new std::mutex[mxOffset + 1];
  std::thread *th = new std::thread[threads - 1];
  int index = 0;
  for (int i = 0; i < threads - 1; ++i) {
    th[i] = std::thread(&SIEFormPMCHW::FillLSEMatrixParallel, this, lseMatrix,
                        std::ref(rwgPerTri), std::ref(index), mut);
  }
  FillLSEMatrixParallel(lseMatrix, rwgPerTri, index, mut);
  for (int i = 0; i < threads - 1; ++i) {
    th[i].join();
  }
  delete[] th;
  delete[] mut;

  std::cout << std::endl;
#ifdef DEBUG
  for (int i = 0; i < 2 * mxOffset; ++i) {
    for (int j = 0; j < 2 * mxOffset; ++j) {
      std::cout << "M(" << i << "," << j << ") = " << std::setprecision(16)
                << (*lseMatrix)(i, j) << "\n";
    }
    std::cout << "\n";
  }
#endif
  return 0;
}

// Worker: fill blocks for one triangle-group i against all j.
int SIEFormPMCHW::FillLSEMatrixParallel(
    blitz::Array<dcmplx, 2> *lseMatrix,
    const std::vector<std::vector<RWGFun *>> &rwgPerTri, int &index,
    std::mutex *mut) {
  using std::vector;
  using std::pair;
  int mxOffset(lseMatrix->shape()(0) / 2);
  while (1) {
    int i;
    {
      // Important : do not remove braces, they are there to block-scope mutex
      // lock
      std::lock_guard<std::mutex> lock(mut[mxOffset]);
      i = index;
      if (i >= rwgPerTri.size()) break;
      ++index;
    }
    std::cout << ".";
    std::cout.flush();
    const auto &trii = rwgPerTri[i];
    if (grnFunLayered == nullptr) // Homogeneous medium
    {
      dcmplx factor1 = I * omega * MU0 * mu;
      dcmplx factor2 = -I * omega * EPS0 * epsilon;
      for (const auto &trij : rwgPerTri) {
        vector<vector<pair<dcmplx, dcmplx>>> DK(
            trii.size(), vector<pair<dcmplx, dcmplx>>(
                            trij.size(), make_pair(dcmplx(0), dcmplx(0))));
        grnFun->SameTriDK(trii, trij, DK);

        for (int k = 0; k < (int)rwgPerTri[i].size(); ++k) {
          int m(trii[k]->EdgeIndex());
          std::lock_guard<std::mutex> lock(mut[m]);
          for (int l = 0; l < (int)trij.size(); ++l) {
            int n(trij[l]->EdgeIndex());
            (*lseMatrix)(m, n) += factor1 * DK[k][l].first;
            (*lseMatrix)(m, n + mxOffset) += DK[k][l].second;
            (*lseMatrix)(m + mxOffset, n) += DK[k][l].second;
            (*lseMatrix)(m + mxOffset, n + mxOffset) += factor2 * DK[k][l].first;
          }
        }
      }
    } else { // Layered medium
#ifdef DEBUG
      vector<vector<pair<dcmplx, dcmplx>>> DK_test(
            trii.size(), vector<pair<dcmplx, dcmplx>>(
                            rwgPerTri[0].size(), make_pair(dcmplx(0), dcmplx(0))));
      grnFunLayered->SameTriDK(trii, rwgPerTri[0], DK_test);
#endif
      for (const auto &trij : rwgPerTri) {
        vector<vector<pair<dcmplx, dcmplx>>> DK_hom(
            trii.size(), vector<pair<dcmplx, dcmplx>>(
                            trij.size(), make_pair(dcmplx(0.0), dcmplx(0.0))));
        
        vector<vector<pair<dcmplx, dcmplx>>> DK_E(
            trii.size(), vector<pair<dcmplx, dcmplx>>(
                            trij.size(), make_pair(dcmplx(0.0), dcmplx(0.0))));

        vector<vector<pair<dcmplx, dcmplx>>> DK_H(
            trii.size(), vector<pair<dcmplx, dcmplx>>(
                            trij.size(), make_pair(dcmplx(1.0), dcmplx(0.0))));
        
        rvec iCenter = trii[0]->TrianglePtr()->Center();
        rvec jCenter = trij[0]->TrianglePtr()->Center();
        if (GetLayerIndex(iCenter) == GetLayerIndex(jCenter)) {
          grnFun[GetLayerIndex(iCenter)].SameTriDK(trii, trij, DK_hom);
        }

        // Calculate the factors based on the layered medium properties
        dcmplx factor1 = I * omega * MU0 * muLayered[GetLayerIndex(jCenter)];
        dcmplx factor2 = -I * omega * EPS0 * epsilonLayered[GetLayerIndex(jCenter)];

        // Calculate DK_E and DK_H for the layered medium
        grnFunLayered->SameTriDK(trii, trij, DK_E);
        grnFunLayered->SameTriDK(trii, trij, DK_H);

        for (int k = 0; k < (int)rwgPerTri[i].size(); ++k) {
          int m(trii[k]->EdgeIndex());
          std::lock_guard<std::mutex> lock(mut[m]);
          for (int l = 0; l < (int)trij.size(); ++l) {
            int n(trij[l]->EdgeIndex());
            (*lseMatrix)(m, n) += factor1 * (DK_hom[k][l].first + DK_E[k][l].first);
            (*lseMatrix)(m, n + mxOffset) += DK_hom[k][l].second + DK_E[k][l].second;
            (*lseMatrix)(m + mxOffset, n) += DK_hom[k][l].second + DK_H[k][l].second;
            (*lseMatrix)(m + mxOffset, n + mxOffset) += 
                                  factor2 * (DK_hom[k][l].first + DK_H[k][l].first);
          }
        }
      }
    }    
  }
  return 0;
}

// ---------------- Matrix + gradient (homogeneous) ----------------
/**
 * @details Triangles are grouped, and workers call FillLSEMatrixAndGradientParallel.
 * Gradient contributions are computed only for boundary-type pairs that
 * cross the 0–1 interface (as per existing code logic).
 */
int SIEFormPMCHW::FillLSEMatrixAndGradient(
    blitz::Array<dcmplx, 2> *lseMatrix,
    blitz::Array<dcmplx, 2> *lseGradientMatrix[3]) {
  using std::vector;
  using std::pair;
  vector<RWGFun>::iterator rwgM, rwgN;

  // Partition the RWG functions to groups over same triangle
  vector<vector<RWGFun *>> rwgPerTri;
  // Whether each triangle is between
  // 1 : 0 and 1
  // 2 : 0 and some other domain
  // 0 : any other case
  vector<int> bordertype;
  for (rwgM = rwgFuns->begin(); rwgM != rwgFuns->end();) {
    rwgPerTri.push_back(vector<RWGFun *>());
    if (std::min(rwgM->TrianglePtr()->BorderingIndeces().first,
                 rwgM->TrianglePtr()->BorderingIndeces().second) == 0) {
      if (std::max(rwgM->TrianglePtr()->BorderingIndeces().first,
                   rwgM->TrianglePtr()->BorderingIndeces().second) == 1) {
        bordertype.push_back(1);
      } else {
        bordertype.push_back(2);
      }
    } else {
      bordertype.push_back(0);
    }
    for (rwgN = rwgM;
         rwgN != rwgFuns->end() and rwgN->TrianglePtr() == rwgM->TrianglePtr();
         rwgN++) {
      rwgPerTri.back().push_back(&(*rwgN));
    }
    rwgM = rwgN;
  }

#ifdef DEBUG
  for (int i = 0; i < rwgTriLen; ++i) {
    std::cout << "Group " << i << " : size " << rwgPerTri[i].size() << "\n";
  }
#endif

  std::cout << "Filling matrix with " << rwgFuns->size() << " elements."
            << std::endl;

  int mxOffset(lseMatrix->shape()(0) / 2);
  std::mutex *mut = new std::mutex[mxOffset + 1];
  std::thread *th = new std::thread[threads - 1];
  int index = 0;
  for (int i = 0; i < threads - 1; ++i) {
    th[i] = std::thread(&SIEFormPMCHW::FillLSEMatrixAndGradientParallel, this,
                        lseMatrix, lseGradientMatrix, std::ref(rwgPerTri),
                        std::ref(bordertype), std::ref(index), mut);
  }
  FillLSEMatrixAndGradientParallel(lseMatrix, lseGradientMatrix, rwgPerTri,
                                   bordertype, index, mut);
  for (int i = 0; i < threads - 1; ++i) {
    th[i].join();
  }
  delete[] th;
  delete[] mut;

  std::cout << std::endl;
#ifdef DEBUG
  for (int i = 0; i < 2 * mxOffset; ++i) {
    for (int j = 0; j < 2 * mxOffset; ++j) {
      std::cout << "M(" << i << "," << j << ") = " << (*lseMatrix)(i, j)
                << "\n";
      for (int k = 0; k < 3; ++k) {
        std::cout << "GM" << k << "(" << i << "," << j
                  << ") = " << (*(lseGradientMatrix[k]))(i, j) << "\n";
      }
    }
    std::cout << "\n";
  }
#endif
  return 0;
}

/**
 * @details Gradient terms are only accumulated for triangle-pair types that meet the
 * selector (bordertype[i] + bordertype[j] == 3), following existing logic.
 */
int SIEFormPMCHW::FillLSEMatrixAndGradientParallel(
    blitz::Array<dcmplx, 2> *lseMatrix,
    blitz::Array<dcmplx, 2> *lseGradientMatrix[3],
    const std::vector<std::vector<RWGFun *>> &rwgPerTri,
    const std::vector<int> &bordertype, int &index, std::mutex *mut) {
  using std::vector;
  using std::pair;
  dcmplx factor1 = I * omega * MU0 * mu;
  dcmplx factor2 = -I * omega * EPS0 * epsilon;
  int mxOffset(lseMatrix->shape()(0) / 2);
  while (1) {
    int i;
    {
      // Important : do not remove braces, they are there to block-scope mutex
      // lock
      std::lock_guard<std::mutex> lock(mut[mxOffset]);
      i = index;
      if (i >= rwgPerTri.size()) break;
      ++index;
    }
    std::cout << ".";
    std::cout.flush();
    const auto &trii = rwgPerTri[i];
    for (int j = 0; j < rwgPerTri.size(); ++j) {
      const auto &trij = rwgPerTri[j];
      // Gradient calculation required only when one triangle is between 0-1 and
      // the other is between 0-!1
      if (bordertype[i] + bordertype[j] != 3) {
        vector<vector<pair<dcmplx, dcmplx>>> DK(
            trii.size(), vector<pair<dcmplx, dcmplx>>(
                             trij.size(), make_pair(dcmplx(0), dcmplx(0))));
        grnFun->SameTriDK(trii, trij, DK);

        for (int k = 0; k < (int)trii.size(); ++k) {
          int m(trii[k]->EdgeIndex());
          std::lock_guard<std::mutex> lock(mut[m]);
          for (int l = 0; l < (int)trij.size(); ++l) {
            int n(trij[l]->EdgeIndex());
            (*lseMatrix)(m, n) += factor1 * DK[k][l].first;
            (*lseMatrix)(m, n + mxOffset) += DK[k][l].second;
            (*lseMatrix)(m + mxOffset, n) += DK[k][l].second;
            (*lseMatrix)(m + mxOffset, n + mxOffset) +=
                factor2 * DK[k][l].first;
          }
        }
      } else {
        double sgn = (bordertype[i] == 1) ? 1 : -1;
        vector<vector<pair<dcmplx, dcmplx>>> DK(
            trii.size(), vector<pair<dcmplx, dcmplx>>(
                             trij.size(), make_pair(dcmplx(0), dcmplx(0))));
        vector<vector<pair<cvec, cvec>>> DKG(
            trii.size(),
            vector<pair<cvec, cvec>>(trij.size(),
                                     make_pair(cvec(0, 0, 0), cvec(0, 0, 0))));
        grnFun->SameTriDKG(trii, trij, DK, DKG);

        for (int k = 0; k < (int)trii.size(); ++k) {
          int m(trii[k]->EdgeIndex());
          std::lock_guard<std::mutex> lock(mut[m]);
          for (int l = 0; l < (int)trij.size(); ++l) {
            int n(trij[l]->EdgeIndex());
            (*lseMatrix)(m, n) += factor1 * DK[k][l].first;
            (*lseMatrix)(m, n + mxOffset) += DK[k][l].second;
            (*lseMatrix)(m + mxOffset, n) += DK[k][l].second;
            (*lseMatrix)(m + mxOffset, n + mxOffset) +=
                factor2 * DK[k][l].first;
            for (int p = 0; p < 3; ++p) {
              (*(lseGradientMatrix[p]))(m, n) +=
                  sgn * factor1 * DKG[k][l].first[p];
              (*(lseGradientMatrix[p]))(m, n + mxOffset) +=
                  sgn * DKG[k][l].second[p];
              (*(lseGradientMatrix[p]))(m + mxOffset, n) +=
                  sgn * DKG[k][l].second[p];
              (*(lseGradientMatrix[p]))(m + mxOffset, n + mxOffset) +=
                  sgn * factor2 * DKG[k][l].first[p];
            }
          }
        }
      }
    }
  }
  return 0;
}

// ----------------------- RHS vector -----------------------
/**
 * @details For homogeneous media, uses IntegrateE/IntegrateH; for layered media,
 * uses IntegrateLayeredEandH to get matching contributions per RWG edge.
 */
int SIEFormPMCHW::FillLSEVector(blitz::Array<dcmplx, 1> *lseVector) {
  std::cout << "Filling vector." << std::endl;
  int mxOffset(lseVector->shape()(0) / 2);
  for (std::vector<RWGFun>::iterator rwgM = rwgFuns->begin();
       rwgM != rwgFuns->end(); rwgM++) {
    for (std::vector<IncidentField *>::iterator incIter = incFields->begin();
         incIter != incFields->end(); incIter++) {
      int m(rwgM->EdgeIndex());
      if (grnFunLayered == nullptr) { // Homogeneous medium
        dcmplx temp;
        temp = (*incIter)->IntegrateE(&(*rwgM));
#ifdef DEBUG
        std::cout << "δ (" << m << ")=" << temp << "\n";
#endif
        (*lseVector)(m) += temp;

        temp = (*incIter)->IntegrateH(&(*rwgM));
#ifdef DEBUG
        std::cout << "κ (" << m << ")=" << temp << "\n";
#endif
        (*lseVector)(m + mxOffset) -= temp;
      } else { // Layered medium
        std::vector<dcmplx> temp;
        temp = (*incIter)->IntegrateLayeredEandH(&(*rwgM));
#ifdef DEBUG
        std::cout << "δ (" << m << ")=" << temp[0] << "\n";
#endif
        (*lseVector)(m) += temp[0];
#ifdef DEBUG
        std::cout << "κ (" << m << ")=" << temp[1] << "\n";
#endif
        (*lseVector)(m + mxOffset) -= temp[1];
      }
    }
#ifdef DEBUG
    std::cout << "\n";
#endif
  }
  return 0;
}

// ---------------- Layer indexing (utility) ----------------
/**
 * @details Indexing convention:
 *  - 0      : below the first interface
 *  - 1..N-1 : between interfaces i-1 and i
 *  - N      : above the last interface
 * with a small tolerance on comparisons.
 */
int SIEFormPMCHW::GetLayerIndex(rvec pos) const {
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

// ---------------- Threading control ----------------
int SIEFormPMCHW::threads = 1;

// Set the number of worker threads used in parallel assembly.
void SIEFormPMCHW::AssignThreads(int t) { threads = t; }
