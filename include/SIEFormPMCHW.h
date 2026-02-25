/**
 * @file SIEFormPMCHW.h
 * @brief Header file for the SIEFormPMCHW class.
 * 
 * PMCHWT (Poggio–Miller–Chang–Harrington–Wu-Tsai) surface-integral formulation
 * for penetrable objects. Supports both homogeneous and layered media.
 *
 * Responsibilities:
 *  - Assemble the linear system matrix A and RHS vector b for PMCHWT.
 *  - Provide parallel assembly helpers (thread-partitioned over triangles).
 *  - In layered media, combine free-space and layered Green's contributions.
 */

// Only include this once during compiling
#ifndef SIEFORMPMCHW_H
#define SIEFORMPMCHW_H

#include <blitz/array.h>
#include <mutex>
#include "GreenF.h"
#include "IncidentField.h"
#include "RWGFun.h"
#include "SIEForm.h"
#include "globals.h"

/// \ingroup formulation
/**
 * @class SIEFormPMCHW
 * @brief PMCHWT formulation for SIE in homogeneous and layered media.
 *
 * Uses RWG basis/testing and Green's function back-ends (homogeneous and
 * optionally layered) to assemble the block matrix:
 *
 * \f[
 * \begin{bmatrix}
 * \mathrm{j}\,\omega\,\mu\,\mathbf{D} & \,\mathbf{K} \\
 * -\,\mathbf{K} & \mathrm{j}\,\omega\,\varepsilon\,\mathbf{D}
 * \end{bmatrix}
 * \begin{bmatrix}
 * \mathbf{J} \\[2pt] \mathbf{M}
 * \end{bmatrix}
 * =
 * \begin{bmatrix}
 * \mathbf{q}^\mathbf{E}_{\mathrm{inc}} \\[2pt] \mathbf{q}^\mathbf{H}_{\mathrm{inc}}
 * \end{bmatrix},
 * \f]
 *
 * where blocks are populated by double-surface integrals of Green kernels
 * and their curls/gradients, evaluated with RWG functions on mesh triangles.
 */
class SIEFormPMCHW : public SIEForm {
 protected:
  // ---- Homogeneous medium parameters ----
  dcmplx epsilon;              ///< Relative permittivity for homogeneous case
  dcmplx mu;                   ///< Relative permeability for homogeneous case

  /** 
   * @brief Wave impedance in homogeneous medium.
   * @param e Relative permittivity.
   * @param m Relative permeability.
   * @return Wave impedance.
   */
  dcmplx Zr(dcmplx e, dcmplx m);

  // ---- Layered medium parameters ----
  std::vector<dcmplx> epsilonLayered;   ///< Layer-wise relative permittivities
  std::vector<dcmplx> muLayered;        ///< Layer-wise relative permeabilities
  std::vector<double> zValsInterfaces;  ///< Interface z-positions (ascending)
  int midLayerIndex;                    ///< Layer index of the mesh middle point.
  GreenF* grnFunLayered;                ///< Layered Green's function provider

  /**
   * @brief Layer-wise wave impedances Z_r.
   * @param e Vector of layer permittivities.
   * @param m Vector of layer permeabilities.
   * @return Vector of impedances per layer.
   */
  std::vector<dcmplx> ZrLayered(const std::vector<dcmplx>& e, 
                                const std::vector<dcmplx>& m);

  /// @brief Number of worker threads for parallel assembly.
  static int threads;

 public:
  // ---------------- Constructors ----------------
  /**
   * @brief Construct PMCHWT for a homogeneous medium.
   * @param inRWGFuns       RWG basis functions (owned elsewhere).
   * @param inGrnFun        Homogeneous Green's function.
   * @param inIncFields     Incident fields used to build the RHS.
   * @param inVacWavelength Vacuum wavelength.
   * @param inEpsilon       Relative permittivity.
   * @param inMu            Relative permeability.
   */
  SIEFormPMCHW(std::vector<RWGFun> *inRWGFuns, GreenF *inGrnFun,
               std::vector<IncidentField *> *inIncFields,
               dcmplx inVacWavelength, dcmplx inEpsilon, dcmplx inMu);
  
  /**
   * @brief Construct PMCHWT for a layered medium.
   * @param inRWGFuns       RWG basis functions (owned elsewhere).
   * @param inGrnFun        Homogeneous Green's function (for same-layer cases).
   * @param inIncFields     Incident fields used to build the RHS.
   * @param inVacWavelength Vacuum wavelength.
   * @param inGrnFunLayered Layered Green's function back-end.
   * @param inEpsilon       Layer-wise relative permittivity values.
   * @param inMu            Layer-wise relative permeability values.
   * @param inZVals         Interface z-positions (ascending).
   * @param inMidLayerIndex Layer index of the mesh middle point.
   */
  SIEFormPMCHW(std::vector<RWGFun> *inRWGFuns,
               GreenF *inGrnFun,
               std::vector<IncidentField *> *inIncFields,
               dcmplx inVacWavelength,
               GreenF *inGrnFunLayered,
               const std::vector<dcmplx>& inEpsilon,
               const std::vector<dcmplx>& inMu,
               const std::vector<double>& inZVals,
               int inMidLayerIndex = 0);               

  // --------------------- Matrix assembly ---------------------
  /**
   * @brief Assemble the PMCHWT system matrix for a homogeneous/layered medium.
   * @param lseMatrix Output A (2Nx2N). Caller must size appropriately.
   * @return 0 on success.
   */
  int FillLSEMatrix(blitz::Array<dcmplx, 2> *lseMatrix);

  /**
   * @brief Worker for parallel homogeneous/layered media assembly (triangle-group partition).
   * @param lseMatrix Output A (shared across workers; protected by row mutexes).
   * @param rwgPerTri RWG basis grouped by owner triangle.
   * @param index     Shared index into @p rwgPerTri (atomic via mutex).
   * @param mut       Array of mutexes; last one guards @p index.
   * @return 0 on success.
   */
  int FillLSEMatrixParallel(blitz::Array<dcmplx, 2> *lseMatrix,
                            const std::vector<std::vector<RWGFun *>> &rwgPerTri,
                            int &index, std::mutex *mut);

  /**
   * @brief Assemble A and its gradients w.r.t. 3 geometry params (homogeneous).
   * @param lseMatrix Output A.
   * @param lseGradientMatrix Array of 3 gradient matrices.
   * @return 0 on success.
   */
  int FillLSEMatrixAndGradient(blitz::Array<dcmplx, 2> *lseMatrix,
                               blitz::Array<dcmplx, 2> *lseGradientMatrix[3]);

  /**
   * @brief Worker for parallel assembly of A and gradients (homogeneous).
   * @param lseMatrix Output A (shared).
   * @param lseGradientMatrix Output gradient matrices (shared).
   * @param rwgPerTri Grouped RWG per triangle.
   * @param bordertype Per-triangle boundary category used for gradient logic.
   * @param index Shared index into @p rwgPerTri (guarded).
   * @param mut Row mutex array; last element guards @p index.
   * @return 0 on success.
   */
  int FillLSEMatrixAndGradientParallel(
      blitz::Array<dcmplx, 2> *lseMatrix,
      blitz::Array<dcmplx, 2> *lseGradientMatrix[3],
      const std::vector<std::vector<RWGFun *>> &rwgPerTri,
      const std::vector<int> &bordertype, int &index, std::mutex *mut);

  // --------------------- Vector assembly ---------------------
  /**
   * @brief Assemble PMCHWT right-hand side vector for current incident fields.
   * @param lseVector Output b (size 2Nx1).
   * @return 0 on success.
   */
  int FillLSEVector(blitz::Array<dcmplx, 1> *lseVector);

  // ------------------------ Utilities ------------------------
  /**
   * @brief Return the layer index for a 3D position (based on z).
   * @param pos Cartesian position.
   * @return Layer index in [0...N], using @ref zValsInterfaces with tolerance.
   */
  int GetLayerIndex(rvec pos) const;
  
  /**
   * @brief Configure how many threads to use for parallel assembly.
   * @param t Number of threads (>=1).
   */
  static void AssignThreads(int t);
};

#endif
