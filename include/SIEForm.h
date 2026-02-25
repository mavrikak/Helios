/**
 * @file SIEForm.h
 * @brief Header for the SIEForm class. The SIEForm class is an abstract base 
 * for Surface Integral Equation (SIE) formulations. It holds shared context 
 * (RWG functions, Green's function, incident fields, angular frequency) and 
 * defines the interface that concrete formulations must implement to assemble 
 * the linear system.
 */

// Only include this once during compiling
#ifndef SIEFORM_H
#define SIEFORM_H

#include <blitz/array.h>
#include "globals.h"

#include <vector>
#include "GreenF.h"
#include "IncidentField.h"
#include "RWGFun.h"

/// \ingroup formulation
/**
 * @class SIEForm
 * @brief Abstract base class for SIE formulations.
 *
 * Derived classes (e.g., PMCHWT) implement the routines that
 * assemble the system matrix and right-hand side using the provided RWG basis,
 * Green's function back-end, and a list of incident fields.
 */
class SIEForm {
 protected:
  std::vector<RWGFun> *rwgFuns;             ///< Pointer to collection of RWG basis functions
  GreenF *grnFun;                           ///< Green's function provider (homogeneous/layered)
  std::vector<IncidentField *> *incFields;  ///< Incident field(s) used to build excitation vector
  dcmplx omega;                             ///< Angular frequency

 public:
  /**
   * @brief Construct the SIE formulation context.
   *
   * Stores pointers to basis functions, Green's function back-end, and
   * incident fields. Computes the angular frequency from the given vacuum
   * wavelength.
   *
   * @param inRWGFuns       RWG basis function array (owned elsewhere).
   * @param inGrnFun        Green's function object (owned elsewhere).
   * @param inIncFields     List of incident fields (owned elsewhere).
   * @param inVacWavelength Vacuum wavelength used to compute omega.
   */
  SIEForm(std::vector<RWGFun> *inRWGFuns, GreenF *inGrnFun,
          std::vector<IncidentField *> *inIncFields, dcmplx inVacWavelength);
  
  /// @brief Virtual destructor.
  virtual ~SIEForm();

  /**
   * @brief Assemble the system matrix @f$ \mathbf{A} @f$.
   *
   * Populates the complex matrix used in the linear system 
   * @f$ \mathbf{A}\mathbf{x}=\mathbf{b} @f$.
   * Derived classes must resize/fill @p matrix entirely.
   *
   * @param matrix Output pointer to a 2D Blitz array to be filled.
   * @return 0 on success; non-zero on error.
   */
  virtual int FillLSEMatrix(blitz::Array<dcmplx, 2> *matrix) = 0;

  /**
   * @brief Assemble system matrix and its gradient w.r.t. geometry parameters.
   *
   * Some inverse/design workflows require gradients of the operator with
   * respect to three geometry parameters (e.g., dx, dy, dz). Derived classes
   * must fill @p matrix and the three gradient matrices in @p gradmatrix.
   *
   * @param matrix     Output system matrix.
   * @param gradmatrix Array of three pointers to gradient matrices 
   *                   @f$(\partial Z / \partial p_i)@f$.
   * @return 0 on success; non-zero on error.
   */
  virtual int FillLSEMatrixAndGradient(blitz::Array<dcmplx, 2> *matrix,
                                       blitz::Array<dcmplx, 2> *gradmatrix[3]) = 0;

  /**
   * @brief Assemble the right-hand side vector @f$ \mathbf{b} @f$.
   *
   * Uses the incident fields and testing with the RWG basis to construct
   * the excitation vector compatible with the chosen formulation.
   *
   * @param vector Output pointer to a 1D Blitz array to be filled.
   * @return 0 on success; non-zero on error.
   */                                      
  virtual int FillLSEVector(blitz::Array<dcmplx, 1> *vector) = 0;
};

#endif
