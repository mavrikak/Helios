/**
 * @file DomainLayered3D.h
 * @brief Layered 3D simulation domain (derived from Domain).
 *
 * A Domain that supports a stack of planar layers with piecewise‑homogeneous,
 * isotropic materials. It owns an array of homogeneous Green's functions (one
 * per layer) used for intra‑layer  interactions and a layered Green's function 
 * used to account for reflections/transmissions across the stack (Sommerfeld 
 * integrals with tabulation). Utilities in #LayeredMediaUtils provide Fresnel 
 * coefficients, transfer matrices, and other per‑layer helpers reused here.
 */

#ifndef DOMAINLAYERED3D_H
#define DOMAINLAYERED3D_H

#include "Domain.h"
#include "GreenFHom3D.h"
#include "GreenFLayered3D.h"
#include "SurfaceMesh.h"
#include <blitz/array.h>
#include <algorithm>
#include <unordered_set>

/// \ingroup domain
/**
 * @class DomainLayered3D
 * @brief 3D stratified (layered) medium domain.
 *
 * Responsibilities on top of #Domain:
 * - Keep per‑layer material parameters and interface z‑positions.
 * - Build per‑layer homogeneous Green's functions plus a layered Green's function
 * with Sommerfeld integration (tabulated on problem‑dependent grids).
 * - Partition arbitrary point sets into layer‑consistent groups (slice).
 * - Compute tabulation ranges and grids (calculateRange, computeGrid).
 * - Provide layered incident sources (PlaneWave, Dipole, Gaussian).
 * - Post‑processing: scattered/total fields and cross sections with layered Green's function.
 */
class DomainLayered3D : public Domain {
 protected:
  // --- Layer data -------------------------------------------------------------
  std::vector<dcmplx> k_L;                  //!< Wave vectors for different layers.
  std::vector<dcmplx> epsilon;              //!< Relative permittivity values for each layer.
  std::vector<dcmplx> mu;                   //!< Relative permeability values for each layer.
  std::vector<double> zValsInterfaces;      //!< Z-coordinates of interface layers.
  std::vector<double> thickness;            //!< Thickness of each internal layer.
  int midLayerIndex;                        //!< Layer index of the mesh middle point.

  // --- Tabulation grids -------------------------------------------------------
  std::vector<Grid> tabulationGrids;        //!< Grid points for tabulation.
  std::vector<Grid> tabulationGridsFields;  //!< Grid points for tabulation (post-processing).

  // --- Green's functions and helpers ------------------------------------------
  GreenF *grnFunLayered;                    //!< Pointer to Green's function for layered media.
  GreenF *grnFunLayeredFields;              //!< Pointer to Green's function for layered media fields calculations.
  LayeredMediaUtils* layeredUtils;          //!< Utility functions for layered media calculations.

  // --- Geometry bookkeeping ---------------------------------------------------
  std::vector<rvec> triPoints;              //!< Points from triangles used for tabulation.
  bool enableLineInt;                       //!< Variable to enable line integrals.

  /**
   * @brief Initialize per‑layer homogeneous Green's functions and the layered Green's function.
   * @return 0 on success.
   */
  virtual int InitGrnFun();



 public:
  // ---------------------------------------------------------------------------
  // Construction & lifetime
  // ---------------------------------------------------------------------------
  /**
   * @brief Constructs a 3D layered domain with given material properties and geometry.
   *
   * This constructor initializes the domain with the provided surface mesh, domain index,
   * permittivity, permeability, and interface positions. It also calculates the layer thicknesses,
   * reads positional data, slices points, computes the tabulation grids, and initializes the Green's function.
   *
   * @param inMesh Pointer to the surface mesh.
   * @param inDIndex Domain index.
   * @param inEpsilon Vector of relative permittivity values for each layer.
   * @param inMu Vector of relative permeability values for each layer.
   * @param inZVals Vector of Z-coordinates defining the layer interfaces.
   * @param inVacuumWavelength The wavelength in vacuum.
   */
  DomainLayered3D(SurfaceMesh *inMesh, int inDIndex, 
                  const std::vector<dcmplx>& inEpsilon, 
                  const std::vector<dcmplx>& inMu, 
                  const std::vector<double>& inZVals, 
                  dcmplx inVacuumWavelength);

  /**
   * @brief Construct without creating Green's functions (owner/aggregate path).
   * @param inMesh Pointer to the surface mesh.
   * @param inDIndex Domain index.
   * @param inEpsilon Vector of relative permittivity values for each layer.
   * @param inMu Vector of relative permeability values for each layer.
   * @param inZVals Vector of Z-coordinates defining the layer interfaces.
   * @param inVacuumWavelength The wavelength in vacuum.
   * @param parent If true, skip InitGrnFun(); caller is responsible to call it later.
   */
  DomainLayered3D(SurfaceMesh *inMesh, int inDIndex, 
                  const std::vector<dcmplx>& inEpsilon, 
                  const std::vector<dcmplx>& inMu, 
                  const std::vector<double>& inZVals, 
                  dcmplx inVacuumWavelength, bool parent);
  
  /** @brief Destructor: deletes layered/homogeneous GreenF objects and utilities. */
  ~DomainLayered3D();

  // ---------------------------------------------------------------------------
  // Domain traits & formulation
  // ---------------------------------------------------------------------------
  /** @brief Always true for this class.
   *  @return true. 
   */
  virtual bool IsLayered() { return true; }

  /**
  * @brief Initialize the layered Green's function used for field evaluation.
  * @param posvec Cloud of observation points; only those inside the domain are used.
  * @details Builds a new set of tabulation grids adapted to @p posvec, creates
  * a dedicated layered Green's function for field evaluation, and wires any dipole
  * incidents to their own post‑processing Green's function.
  */
  virtual void newGrnFunLayered(std::vector<rvec> &posvec);

  /**
   * @brief Initialize the SIE formulation (PMCHWT) for the given wavelength.
   * @param inVacWavelength Vacuum wavelength \f$\lambda_0\f$ used to configure the formulation.
   * @return 0 on success; non-zero otherwise.
   */
  int InitFormulation(dcmplx inVacWavelength);

  // ---------------------------------------------------------------------------
  // Layer‑optics helpers (Fresnel/transfer, secondary fields, 1D tests)
  // ---------------------------------------------------------------------------
  /**
   * @brief Computes Fresnel reflection and transmission coefficients for an interface.
   * @param layerIndex1 Index of the first layer.
   * @param layerIndex2 Index of the second layer.
   * @param krho Transverse wave number.
   * @param waves Wave type (TE or TM).
   * @return Vector of complex coefficients.
   */
  std::vector<dcmplx> FresnelCoeff(int layerIndex1, int layerIndex2, 
                                   dcmplx krho, const std::string& waves) const;

  /**
   * @brief Computes the transfer matrix for an interface.
   * @param layerIndexI First layer index.
   * @param layerIndexJ Second layer index.
   * @param krho Transverse wave number.
   * @param waves Wave type (TE or TM).
   * @return Transfer matrix as a 2D complex vector.
   */
  std::vector<std::vector<dcmplx>> InterfaceTransferMatrix(
    int layerIndexI, int layerIndexJ, dcmplx krho, const std::string& waves) const;

  /**
   * @brief Computes secondary wave coefficients.
   * @param krho Transverse wave number.
   * @param observationLayer Layer of observation.
   * @param sourceLayer Layer of source.
   * @param waves Wave type (TE or TM).
   * @return 2D array of complex coefficients.
   */
  blitz::Array<dcmplx, 2> Secondary(dcmplx krho, int observationLayer, 
                                    int sourceLayer, const std::string& waves);

  /**
   * @brief Computes fields at observation points.
   * @param krho Transverse wave number.
   * @param observationPointsZ Vector of Z-coordinates for observation points.
   * @param sourcePointZ Z-coordinate of the source point.
   * @param waves Wave type (TE or TM).
   * @return Vector of computed field values.
   */
  std::vector<dcmplx> Fields(dcmplx krho, std::vector<double> observationPointsZ, 
                             double sourcePointZ, const std::string& waves);
  
  /**
   * @brief Compute fields over a 1D sweep and write them to a text file.
   * @param krho Transverse wavenumber.
   * @param start Start value of the sweep parameter.
   * @param end End value of the sweep parameter.
   * @param numPoints Number of uniformly spaced points in the sweep.
   * @param sourcePointZ Source point z-coordinate.
   * @param waves Wave type selector ("TE" or "TM").
   * @return 0 on success; non-zero on I/O or computation error.
   */
  int writeFields(dcmplx krho, double start, double end, int numPoints, 
                  double sourcePointZ, const std::string& waves);

  // ---------------------------------------------------------------------------
  // Material queries and layer indexing
  // ---------------------------------------------------------------------------
  /**
   * @brief Returns the layer index corresponding to a spatial position.
   * @param pos Position in space.
   * @return Layer index.
   */
  int GetLayerIndex(rvec pos) const;

  /**
   * @brief Relative permittivity at a spatial location.
   * @param pos Query location.
   * @return Complex relative permittivity \f$\varepsilon_r(\mathbf r)\f$.
   */
  virtual dcmplx Epsilon(rvec pos);

  /**
   * @brief Relative permeability at a spatial location.
   * @param pos Query location.
   * @return Complex relative permeability \f$\mu_r(\mathbf r)\f$.
   */
  virtual dcmplx Mu(rvec pos);

  /**
   * @brief Accessor to the Green's function used for layered-media post-processing.
   * @return Pointer to the layered Green's function instance.
   */
  virtual GreenF* GetGrnFunLayeredFields() const { return grnFunLayeredFields; }

  /**
   * @brief Access the full vector of layer permittivities.
   * @return Complex relative permittivities \f$\mathbf{\varepsilon_r}(\mathbf r)\f$.
   */
  std::vector<dcmplx> GetEpsilon() const { return epsilon; }

  // ---------------------------------------------------------------------------
  // Position partitioning & tabulation grid generation
  // ---------------------------------------------------------------------------
  /**
   * @brief Computes the slicing of points with respect to interfaces.
   * @param pos1 Primary set of positions.
   * @param pos2_opt Optional secondary set of positions.
   * @return Vector of grouped positions.
   */
  std::vector<PositionGroup> slice(const std::vector<rvec>& pos1,
                                   std::vector<rvec>* pos2_opt = nullptr);
  
  /**
   * @brief Calculates the tabulation range for a group of points.
   * @param pts Group of positions.
   * @param rmin Minimum radial range.
   * @param zmin Minimum vertical range.
   * @param margin Additional margin buffer.
   * @return Computed range object.
   */
  Range calculateRange(const PositionGroup& pts, double rmin = 0.01,
                       double zmin = 0.01, double margin = 0.05);
  
  /**
   * @brief Computes the grid points for tabulation.
   * @param pos Vector of grouped positions.
   * @param nr Number of radial grid points.
   * @param nz Number of vertical grid points.
   * @param rmin Minimum radial range.
   * @param margin Additional margin buffer.
   * @param zmin Minimum vertical range.
   * @return Computed grid points.
   */
  std::vector<Grid> computeGrid(const std::vector<PositionGroup>& pos,
                                int nr = 20, int nz = 20, double rmin = 1e-2,
                                double margin = 0.05, double zmin = 1e-2);
  
  /**
   * @brief Self-test: slicing of points with respect to interfaces (prints to console).
   */
  void test_slice();

  /**
   * @brief Self-test: tabulation range computation (prints to console).
   */
  void test_range();

  /**
   * @brief Self-test: grid generation (prints to console).
   */
  void test_grid();

  // ---------------------------------------------------------------------------
  // Incident fields (layered variants)
  // ---------------------------------------------------------------------------
  /**
   * @brief Add a PlaneWave excitation adapted to the layered medium.
   * @param wavelength Vacuum wavelength \f$\lambda_0\f$ (complex allowed).
   * @param propagationDirection Unit vector \f$\hat{\mathbf k}\f$ (direction of travel).
   * @param polarization Complex polarization vector (transverse to \f$\hat{\mathbf k}\f$).
   * @return 0 on success; non-zero otherwise.
   */
  virtual int AddPlaneWave(dcmplx wavelength, rvec propagationDirection,
                           cvec polarization);

  /**
   * @brief Add a Dipole excitation in the layered medium.
   * @param position Dipole position.
   * @param polarization Complex dipole moment vector.
   * @return 0 on success; non-zero otherwise.
   */
  virtual int AddDipole(rvec position, cvec polarization);

  /**
   * @brief Add a Gaussian beam excitation in the layered medium.
   * @param wavelength Vacuum wavelength.
   * @param waist Beam waist radius \f$ w_0 \f$.
   * @param focal_point Beam focus position.
   * @param propagationdirection Unit vector of beam propagation.
   * @param polarization Complex polarization vector.
   * @return 0 on success; non-zero otherwise.
   * @note Not implemented for layered media.
   */
  virtual int AddGaussian(double wavelength, double waist, rvec focal_point,
                          rvec propagationdirection, cvec polarization);
  
  // ---------------------------------------------------------------------------
  // Field evaluation (uses layered GFs)
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
  virtual std::pair<cvec, cvec> SecFieldEandH(blitz::Array<dcmplx, 1> *solVector, rvec pos);
  
  // ---------------------------------------------------------------------------
  // Cross sections
  // ---------------------------------------------------------------------------
  /**
   * @brief Scattering cross section.
   * @param solVector Solution coefficients \f$\mathbf{x}\f$.
   * @return Scattering cross section and near-field intensity.
   */
  virtual std::vector<double> scatteringCS(blitz::Array<dcmplx, 1> *solVector);

  /**
   * @brief Absorption cross section.
   * @param solVector Solution coefficients \f$\mathbf{x}\f$.
   * @return Absorption cross section.
   * @note Not always accurate (volumetric methods based on Ohmic losses have better accuracy).
   */
  virtual double absorptionCS(blitz::Array<dcmplx, 1> *solVector);

  /**
   * @brief Extinction cross section.
   * @param solVector Solution coefficients \f$\mathbf{x}\f$.
   * @return Extinction cross section.
   * @note Not always accurate (volumetric methods based on Ohmic losses have better accuracy).
   */
  virtual double extinctionCS(blitz::Array<dcmplx, 1> *solVector);
  
  // ---------------------------------------------------------------------------
  // Debug / logging utilities
  // ---------------------------------------------------------------------------
  /**
   * @brief Append a single vector field sample to debug text files.
   * @param fieldName Logical field name (e.g., "E_sc" or "H_inc") used in filenames.
   * @param vec 3-component complex vector sample to append.
   */
  void dumpField(const std::string& fieldName, 
                 const blitz::TinyVector<std::complex<double>, 3>& vec);

  // ---------------------------------------------------------------------------
  // Incident cleanup
  // ---------------------------------------------------------------------------
  /**
   * @brief Clears all incident fields and resets the cached layered Green's 
   * function for field evaluation.
   */
  virtual int ClearIncidentFields() override;
};

#endif
