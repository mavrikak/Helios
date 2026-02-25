/**
 * @file Domain.h
 * @brief Abstract simulation domain class for SIE scattering.
 *
 * A Domain represents a volumetric region bounded by a surface mesh on which
 * Rao–Wilton–Glisson (RWG) basis functions live. It stores incident fields,
 * constructs a BSP tree for robust inside/outside queries, and provides a
 * Green's function handle and a chosen SIE formulation to assemble the linear
 * system. Concrete subclasses (e.g., homogeneous vs. layered) supply material
 * properties and cross‑section computations.
 */


#ifndef DOMAIN_H
#define DOMAIN_H

#include <blitz/array.h>
#include <fstream>
#include "Bsp3D.h"
#include "GreenF.h"
#include "IncidentField.h"
#include "RWGFun.h"
#include "SIEForm.h"
#include "SurfaceMesh.h"
#include "Triangle.h"
#include "globals.h"

/// \ingroup domain
/**
 * @class Domain
 * @brief Abstract base class for domains.
 *
 * Provides methods to assemble the linear system of equations, compute fields,
 * and query material properties. Derived classes must implement Green's
 * function initialization, material permittivity/permeability, and scattering
 * metrics.
 */
class Domain {
 protected:
  // --- Discretization and excitations ---
  std::vector<RWGFun> rwgFuns;                  //!< RWG basis functions on the boundary.
  std::vector<IncidentField *> incidentFields;  //!< List of incident excitations.
  Bsp3D bspTree;                                //!< BSP tree for classification.

  // --- Physical properties and formulations ---
  GreenF *grnFun;                               //!< Green's function for this domain.
  SIEForm *formulation;                         //!< Formulation used for SIE.

  // --- Domain parameters ---
  dcmplx vacuumWavelength;                      //!< Wavelength in vacuum.
  dcmplx omega;                                 //!< Angular frequency.
  int index;                                    //!< Domain index (mesh region).
  static bool NeedsAccurate;                    //!< Flag for accurate computation.
  
  /**
   * @brief Initialize the Green's function (to be implemented by derived classes).
   * @return 0 on success; non-zero otherwise.
   */
  virtual int InitGrnFun() = 0;

 public:
  // ---------------------------------------------------------------------------
  // Construction & lifetime
  // ---------------------------------------------------------------------------
  /**
   * @brief Construct a domain over a mesh region at a given vacuum wavelength.
   * @param mesh Surface mesh containing all boundary triangles.
   * @param dIndex Integer region index (as labeled in @p mesh).
   * @param inVacuumWavelength Vacuum wavelength \f$\lambda_0\f$ (complex to allow dispersion).
   */
  Domain(SurfaceMesh *mesh, int dIndex, dcmplx inVacuumWavelength);

  /** @brief Virtual destructor. */
  virtual ~Domain();

  // ---------------------------------------------------------------------------
  // Mesh / basis assembly
  // ---------------------------------------------------------------------------
   /**
   * @brief Add all boundary triangles belonging to region @p dIndex and build RWGs.
   * @param mesh Input surface mesh.
   * @param dIndex Region index to extract.
   * @return Number of boundary triangles added.
   */
  int AddTrianglesFromMesh(SurfaceMesh *mesh, int dIndex);

  /** @brief Debug: print triangle centroids and RWG evaluations at centroids. 
   * @return 0 on success.
   */
  int PrintAllRWGCentroids();

  /** @brief Debug: print all edges and the RWGs that live on them. 
   * @param mesh Input surface mesh.
   * @return 0 on success.
   */
  int PrintAllEdges(SurfaceMesh *mesh);

  /**
   * @brief Direct access to an RWG basis function by global index.
   * @param index Global RWG index.
   * @return Pointer to the corresponding RWG basis function.
   */
  RWGFun *RGWFunPtr(int index);

  // ---------------------------------------------------------------------------
  // Domain metadata
  // ---------------------------------------------------------------------------
  /**
   * @brief Return this domain's integer index.
   * @return Region index associated with this domain.
   */
  int Index();

  /** @brief Enable globally more accurate (and usually slower) evaluation. */
  static void EnableAccurate();
  
  // ---------------------------------------------------------------------------
  // Formulation assembly hooks
  // ---------------------------------------------------------------------------
  /**
   * @brief Initialize the SIE formulation object at the given wavelength.
   * @param vacWavelength Vacuum wavelength \f$\lambda_0\f$ used to configure the formulation.
   * @return 0 on success; non-zero otherwise.
   */
  virtual int InitFormulation(dcmplx vacWavelength) = 0;

  /**
   * @brief Assemble the linear system matrix using the active formulation.
   * @param matrix Output system matrix \f$\mathbf{A}\f$ (modified in place).
   * @return 0 on success; non-zero otherwise.
   */
  virtual int FillLSEMatrix(blitz::Array<dcmplx, 2> *matrix);

  /**
   * @brief Assemble the linear system matrix and its gradient w.r.t. geometry/parameters.
   * @param matrix Output system matrix \f$\mathbf{A}\f$.
   * @param gradmatrix Array of three pointers to gradient matrices.
   * @return 0 on success; non-zero otherwise.
   */
  virtual int FillLSEMatrixAndGradient(blitz::Array<dcmplx, 2> *matrix,
                                       blitz::Array<dcmplx, 2> *gradmatrix[3]);
                                    
  /**
   * @brief Assemble the right-hand-side (RHS) vector using the active formulation.
   * @param vector Output RHS \f$\mathbf{b}\f$ (modified in place).
   * @return 0 on success; non-zero otherwise.
   */
  virtual int FillLSEVector(blitz::Array<dcmplx, 1> *vector);

  // ---------------------------------------------------------------------------
  // Geometry queries
  // ---------------------------------------------------------------------------
  /**
   * @brief Test whether a spatial point lies inside the domain.
   * @param point Query location.
   * @return true if inside; false otherwise.
   */
  bool IsInside(rvec point);
  
  /** @brief Is this a layered domain? Must be defined by subclasses.
   *  @return true for layered domains; false otherwise. 
   */
  virtual bool IsLayered() = 0;
  
  /**
   * @brief Check whether any incident field is a dipole.
   * @return true if a dipole source exists; false otherwise.
   */
  bool hasDipoles();
  
  // ---------------------------------------------------------------------------
  // Field evaluation (post‑processing)
  // ---------------------------------------------------------------------------
  /**
   * @brief Secondary (scattered) electric field at a point from the solution vector.
   * @param solVector Solution coefficients \f$\mathbf{x}\f$ (current unknowns).
   * @param pos Observation point.
   * @return Complex electric field vector \f$\mathbf{E}^{sc}(\mathbf{r})\f$.
   */
  virtual cvec SecFieldE(blitz::Array<dcmplx, 1> *solVector, rvec pos);
  
  /**
   * @brief Total incident electric field at a point from all incident sources.
   * @param pos Observation point.
   * @return Complex electric field vector \f$\mathbf{E}^{inc}(\mathbf{r})\f$.
   */
  virtual cvec IncFieldE(rvec pos);
  
  /**
   * @brief Secondary (scattered) magnetic field at a point from the solution vector.
   * @param solVector Solution coefficients \f$\mathbf{x}\f$ (current unknowns).
   * @param pos Observation point.
   * @return Complex magnetic field vector \f$\mathbf{H}^{sc}(\mathbf{r})\f$.
   */
  virtual cvec SecFieldH(blitz::Array<dcmplx, 1> *solVector, rvec pos);
  
  /**
   * @brief Total incident magnetic field at a point from all incident sources.
   * @param pos Observation point.
   * @return Complex magnetic field vector \f$\mathbf{H}^{inc}(\mathbf{r})\f$.
   */
  virtual cvec IncFieldH(rvec pos);
  
  /**
   * @brief Secondary fields \f$(\mathbf E^{sc}, \mathbf H^{sc})\f$ at a point.
   * @param solVector Solution coefficients \f$\mathbf{x}\f$.
   * @param pos Observation point.
   * @return Pair \f$\{\mathbf E^{sc}(\mathbf r), \mathbf H^{sc}(\mathbf r)\}\f$.
   */
  virtual std::pair<cvec, cvec> SecFieldEandH(blitz::Array<dcmplx, 1> *solVector,
                                              rvec pos);
  
  // ---------------------------------------------------------------------------
  // Power & cross sections
  // ---------------------------------------------------------------------------
  /**
   * @brief Power crossing the domain boundary computed from the RWG solution.
   * @param solVector Solution coefficients \f$\mathbf{x}\f$.
   * @return Time-averaged outward power.
   * @note Not implemented for layered media.
   */
  double Power(blitz::Array<dcmplx, 1> *solVector);
  
  /**
   * @brief Scattering cross section.
   * @param solVector Solution coefficients \f$\mathbf{x}\f$.
   * @return Scattering cross section and near-field intensity.
   */
  virtual std::vector<double> scatteringCS(blitz::Array<dcmplx, 1> *solVector) = 0;
  
  /**
   * @brief Absorption cross section.
   * @param solVector Solution coefficients \f$\mathbf{x}\f$.
   * @return Absorption cross section.
   * @note Not always accurate (volumetric methods based on Ohmic losses have better accuracy).
   */
  virtual double absorptionCS(blitz::Array<dcmplx, 1> *solVector) = 0;

  /**
   * @brief Extinction cross section.
   * @param solVector Solution coefficients \f$\mathbf{x}\f$.
   * @return Extinction cross section.
   * @note Not always accurate (volumetric methods have better accuracy).
   */
  virtual double extinctionCS(blitz::Array<dcmplx, 1> *solVector) = 0;

  // ---------------------------------------------------------------------------
  // Incident field management
  // ---------------------------------------------------------------------------
  /**
   * @brief Add a new incident source (ownership transferred to Domain).
   * @param inIncident Pointer to dynamically allocated incident field.
   * @return 0 on success; non-zero otherwise.
   */
  int AddIncident(IncidentField *inIncident);

  /** @brief Delete and clear all incident sources.
   *  @return 0 on success. 
   */
  virtual int ClearIncidentFields();

  // ---------------------------------------------------------------------------
  // Material properties (implemented by subclasses)
  // ---------------------------------------------------------------------------
  /**
   * @brief Relative permittivity at a spatial location.
   * @param pos Query location.
   * @return Complex relative permittivity \f$\varepsilon_r(\mathbf r)\f$.
   */
  virtual dcmplx Epsilon(rvec pos) = 0;
  
  /**
   * @brief Relative permeability at a spatial location.
   * @param pos Query location.
   * @return Complex relative permeability \f$\mu_r(\mathbf r)\f$.
   */
  virtual dcmplx Mu(rvec pos) = 0;

  // ---------------------------------------------------------------------------
  // Layered‑media helpers
  // ---------------------------------------------------------------------------
  /**
   * @brief Accessor to the Green's function used for layered-media post-processing.
   * @return Pointer to the layered Green's function instance.
   */
  virtual GreenF* GetGrnFunLayeredFields() const = 0;
  
  /**
   * @brief Build or refresh layered Green's function tables for positions of interest.
   * @param posvec Vector of observation positions used to seed/refresh tables.
   */
  virtual void newGrnFunLayered(std::vector<rvec> &posvec) = 0;

  // ---------------------------------------------------------------------------
  // Convenience incident creators (to be implemented by subclasses)
  // ---------------------------------------------------------------------------
  /**
   * @brief Add a PlaneWave excitation (ownership transferred to Domain).
   * @param wavelength Vacuum wavelength \f$\lambda_0\f$ (complex allowed).
   * @param propagationDirection Unit vector \f$\hat{\mathbf k}\f$ (direction of travel).
   * @param polarization Complex polarization vector (transverse to \f$\hat{\mathbf k}\f$).
   * @return 0 on success; non-zero otherwise.
   */
  virtual int AddPlaneWave(dcmplx wavelength, rvec propagationDirection,
                           cvec polarization) = 0;
  
  /**
   * @brief Add a Dipole excitation (ownership transferred to Domain).
   * @param location Dipole position.
   * @param polarization Complex dipole moment vector.
   * @return 0 on success; non-zero otherwise.
   */
  virtual int AddDipole(rvec location, cvec polarization) = 0;

  /**
   * @brief Add a Gaussian beam excitation.
   * @param wavelength Vacuum wavelength.
   * @param waist Beam waist radius \f$ w_0 \f$.
   * @param focal_point Beam focus position.
   * @param propagationdirection Unit vector of beam propagation.
   * @param polarization Complex polarization vector.
   * @return 0 on success; non-zero otherwise.
   * @note Not implemented for layered media.
   */
  virtual int AddGaussian(double wavelength, double waist, rvec focal_point,
                          rvec propagationdirection, cvec polarization) = 0;
  
  // ---------------------------------------------------------------------------
  // Lookup tables
  // ---------------------------------------------------------------------------
  /**
   * @brief Provide an erfc lookup table to the Green's function implementation.
   * @param tin Pointer to a binary erfc lookup table.
   * @return 0 on success; non-zero otherwise.
   */
  int InitLookupTable(LookupTableBin *tin);
};

#endif
