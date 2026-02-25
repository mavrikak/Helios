/**
* @file GreenFLayered3D.h
* @brief Layered‑medium 3D Green's function (scalar/gradient/dyadic) + RWG integrals.
*
* This header declares a comprehensive Green's function engine for **stratified media**
* (stack of planar layers with piecewise‑constants \f$\epsilon_{r,i},\mu_{r,i}\f$). 
* It supports:
*
* - Pointwise kernels: scalar \f$G\!,\;\nabla'G\!,\;\overline{\overline G}\f$ with
* near‑singular **smoothed** variants for robust area integration.
* - Assembly helpers: two‑triangle (D/K) and point–triangle (delta/kappa) integrals with
* singularity cancellation and optional **edge line integrals** on interfaces.
* - **Tabulation + interpolation** of Sommerfeld integrals on structured grids, with
* per‑interaction coefficient sets (intra‑layer, inter‑layer) to accelerate reuse.
*
* Design notes
* ------------
* - A **SommerfeldIntegrator** is owned per interaction group; it evaluates smooth parts of
* layered Green's terms in cylindrical spectral coordinates. Quasi‑static (QS) pieces are
* handled analytically (cancel singularities), while the radiative remainder is tabulated.
* - The `coeffIntra*`/`coeffInter` vectors hold precomputed interpolation weights to map (r,z)
* to the tabulation grid(s); `intra*Index`/`interIndeces` bind triangles/layers to the
* right coefficient set.
* - `enableLineIntegrals` toggles analytic edge terms used when triangles share an interface.
*/

#ifndef GREENFLAYERED3D_H
#define GREENFLAYERED3D_H

#include <blitz/array.h>
#include <cmath>
#include <tuple>
#include "GaussQuad.h"
#include "GreenF.h"
#include "globals.h"
#include "SommerfeldIntegrator.h"
#include "SingCancellation.h"
#include "TriQuadPol.h"

/// \ingroup greenFunction
/**
 * @class GreenFLayered3D
 * @brief Layered background Green's function with tabulation and D/K/delta/kappa integrals.
 *
 * A high-level facade around Sommerfeld-integral evaluation for layered media.
 * It owns:
 * - Per-layer wavenumbers @p k_L (from GreenF base class).
 * - Medium parameters @p epsilon, @p mu and interface depths
 * @p zValsInterfaces with derived @p thickness.
 * - Precomputed coefficient tables for different interaction types:
 * intra-layer (top/bottom vs internal) and inter-layer.
 * - Interpolation grids (r, z1[, z2]) per interaction group, chosen from
 * mesh-based ranges (triangle barycenters / quadrature points).
 *
 * The class exposes pointwise kernels (scalar/grad/dyadic) and assembly helpers
 * (SameTriDK, SameTriDeltaKappa) used by SIEFormPMCHW.
 */
class GreenFLayered3D : public GreenF {
 private:
  // ------------------------------------------------------------------------------
  // -------------------------- Layer stack & tabulation --------------------------
  // ------------------------------------------------------------------------------
  std::vector<dcmplx> epsilon;                          /*!< Layered medium electric permittivity. */
  std::vector<dcmplx> mu;                               /*!< Layered medium magnetic permeability. */
  std::vector<double> zValsInterfaces;                  /*!< z-coordinates of layered medium interfaces. */
  std::vector<double> thickness;                        /*!< Internal layers' thicknesses. */
  std::vector<Grid> tabulationGrids;                    /*!< Tabulation grids for Green's function. */
  int midLayerIndex;                                    /*!< Layer index of the mesh middle point. */
  LayeredMediaUtils* layeredUtils;                      /*!< Utility functions for layered media domain. */
  
  // Interpolation coeffs (grouped by interaction type)
  std::vector<CoeffsIntra> coeffIntraTB;                /*!< Interpolation coefficients for top/bottom intra-layer interactions. */
  std::vector<CoeffsIntra> coeffIntraMinus;             /*!< Interpolation coefficients for internal intra-layer interactions. */
  std::vector<CoeffsIntra> coeffIntraPlus;              /*!< Interpolation coefficients for internal intra-layer interactions. */
  std::vector<CoeffsInter> coeffInter;                  /*!< Interpolation coefficients for inter-layer interactions. */
  
  // Indices linking triangles/layers to coefficient sets
  std::vector<int> intraTBIndex;                        /*!< Indices for top/bottom intra-layer interactions. */
  std::vector<int> intraPMIndex;                        /*!< Indices for internal intra-layer interactions. */
  std::vector<std::complex<int>> interIndeces;          /*!< Indices for inter-layer interactions. */
  
  // Per‑case Sommerfeld integrators (smooth parts)
  std::vector<SommerfeldIntegrator*> sommIntIntraTB;    /*!< Sommerfeld integrator for top/bottom intra-layer interactions. */
  std::vector<SommerfeldIntegrator*> sommIntIntraMinus; /*!< Sommerfeld integrator for internal intra-layer interactions. */
  std::vector<SommerfeldIntegrator*> sommIntIntraPlus;  /*!< Sommerfeld integrator for internal intra-layer interactions. */
  std::vector<SommerfeldIntegrator*> sommIntInter;      /*!< Sommerfeld integrator for inter-layer interactions. */
  
  bool enableLineIntegrals = false;                     /*!< Flag to enable line integrals. */

  GaussQuad<GreenFLayered3D> gauss;                     /*!< Gaussian quadrature object for integrations. */

 public:
  // -----------------------------------------------------------------------------
  // ------------------------------- Thresholds ----------------------------------
  // -----------------------------------------------------------------------------
  /**
   * @brief Checks whether the given positions exceed a predefined distance threshold.
   * @param r First position.
   * @param rp Second position.
   * @return True if above threshold, false otherwise.
   */
  bool AboveThreshold(rvec r, rvec rp);

  /**
   * @brief Checks whether the given point-triangle pair exceeds a predefined distance threshold.
   * @param r Observation point.
   * @param Tp Source triangle.
   * @return True if above threshold, false otherwise.
   */
  bool AboveThreshold(rvec r, Triangle* Tp);

  /**
   * @brief Checks whether the given triangles exceed a predefined threshold.
   * @param T First triangle.
   * @param Tp Second triangle.
   * @return True if above threshold, false otherwise.
   */
  bool AboveThreshold(Triangle* T, Triangle* Tp);

  // -----------------------------------------------------------------------------
  // ---------------------------- Constructors/destructors -----------------------
  // -----------------------------------------------------------------------------
  /**
   * @brief Construct layered Green's function (non-magnetic layers).
   * @param inWavelength Vacuum wavelength \f$\lambda_0\f$.
   * @param inEpsilonMedium Relative permittivities for all layers (bottom->top).
   * @param inZValsInterfaces Monotone increasing z of interfaces (size = nLayers-1).
   * @param inThickness Layer thicknesses (size = nLayers-2).
   * @param inTabulationGrids Interpolation grids per interaction bucket.
   * @param enableLineInt Enable interface line-integrals for singular cancellation.
   * @param inMidLayerIndex Layer index of the middle point in the mesh.
   */
  GreenFLayered3D(dcmplx inWavelength,
                  const std::vector<dcmplx>& inEpsilonMedium,
                  const std::vector<double>& inZValsInterfaces,
                  const std::vector<double>& inThickness,
                  const std::vector<Grid>& inTabulationGrids,
                  const bool enableLineInt = false,
                  int inMidLayerIndex = 0);
  
  /**
   * @brief Construct layered Green's function (magnetic layers).
   * @param inWavelength Vacuum wavelength \f$\lambda_0\f$.
   * @param inEpsilonMedium Relative permittivities for all layers (bottom -> top).
   * @param inMuMedium Relative permeabilities for all layers (bottom -> top).
   * @param inZValsInterfaces Monotone increasing z of interfaces (size = nLayers-1).
   * @param inThickness Layer thicknesses (size = nLayers-2).
   * @param inTabulationGrids Interpolation grids per interaction bucket.
   * @param enableLineInt Enable interface line-integrals for singular cancellation.
   * @param inMidLayerIndex Layer index of the middle point in the mesh.
   */
  GreenFLayered3D(dcmplx inWavelength,
                  const std::vector<dcmplx>& inEpsilonMedium,
                  const std::vector<dcmplx>& inMuMedium,
                  const std::vector<double>& inZValsInterfaces,
                  const std::vector<double>& inThickness,
                  const std::vector<Grid>& inTabulationGrids,
                  const bool enableLineInt = false,
                  int inMidLayerIndex = 0);
  
  /** @brief Destructor. */
  virtual ~GreenFLayered3D(){ delete layeredUtils; layeredUtils = nullptr; };

  // -----------------------------------------------------------------------------
  // ---------------------------- Geometry helpers -------------------------------
  // -----------------------------------------------------------------------------
  /**
   * @brief Return layer index of a 3D point.
   * @param pos Cartesian position.
   * @return 0 for bottom half-space, 1..n-1 for internal layers, n for top half-space.
   */
  int GetLayerIndex(rvec pos) const;

  /** @brief Writes tabulation grids to a file.
   * @param filename Name of the file.
   * @param grids Tabulation grids.
   */
  void writeGridsToFile(const std::string& filename, const std::vector<Grid>& grids);

  // -----------------------------------------------------------------------------
  // -------------------------- Layered medium helper ----------------------------
  // -----------------------------------------------------------------------------
  /**
   * @brief Getter for layeredUtils pointer.
   * @return Pointer to layeredUtils.
   */
  LayeredMediaUtils* getLayeredUtils() { return layeredUtils; }

  // -----------------------------------------------------------------------------
  // ----------------------------- Pointwise kernels -----------------------------
  // -----------------------------------------------------------------------------
  /**
   * @brief Evaluates the homogeneous Green's function for given positions.
   * @param r Observation position.
   * @param rp Source position.
   * @return Scalar Green's function value.
   */
  virtual dcmplx Evaluate(rvec r, rvec rp);

  /**
   * @brief Evaluates the dyadic homogeneous Green's function.
   * @param r Observation position.
   * @param rp Source position.
   * @return Dyadic Green's function value.
   */
  virtual cdyad EvaluateDyadic(rvec r, rvec rp);

  /**
   * @brief Evaluates the secondary terms for dipole incident field in layered background.
   * @param r Observation position.
   * @param rp Source position.
   * @param dyadicType Type of dyadic Green's function ('normal' or 'curl').
   * @return Dyadic Green's function value.
   */
  cdyad EvaluateDyadic(rvec r, rvec rp, std::string dyadicType);

  /**
   * @brief Computes the derivative with respect to rp (homogeneous).
   * @param r Observation position.
   * @param rp Source position.
   * @return Gradient vector.
   */
  virtual cvec Gradient(rvec r, rvec rp);

  /**
   * @brief Computes all the aforementioned (homogeneous) with respect to rp.
   * @param r Observation position.
   * @param rp Source position.
   * @return Tuple of scalar, gradient vector, dyadic Green's function.
   */
  virtual std::tuple<dcmplx, cvec, cdyad> EvaluateAll(rvec r, rvec rp);

  /**
   * @brief Smoothed homogeneous Green's function.
   * @param r Observation position.
   * @param rp Source position.
   * @return Smoothed scalar Green's function value.
   */
  virtual dcmplx Smoothed(rvec r, rvec rp);

  /**
  * @brief Half-smoothed homogeneous Green's function.
  * @param r Observation position.
  * @param rp Source position.
  * @return Half-smoothed scalar Green's function value.
  */
  virtual dcmplx halfSmoothed(rvec r, rvec rp);
  
  /**
  * @brief Smoothed homogeneous Green's function gradient.
  * @param r Observation position.
  * @param rp Source position.
  * @return Smoothed gradient vector.
  */
  virtual cvec GradientSmoothed(rvec r, rvec rp);

  // -----------------------------------------------------------------------------
  // ------------------------------- Integrals -----------------------------------
  // -----------------------------------------------------------------------------
  /**
   * @brief D-operator matrix element 
   * \f$\langle \mathbf{f}(\mathbf{r}),\,\overline{\overline G},\,\mathbf{f}(\mathbf{r}^{\prime})\rangle\f$.
   * @param f Testing RWG function (on observation triangle).
   * @param fp Basis RWG function (on source triangle).
   * @return Complex scalar entry of \f$\mathbf D\f$.
   * @note Not used here.
   */
  virtual dcmplx IntegrateD(RWGFun* f, RWGFun* fp);
  
  /**
   * @brief K-operator matrix element 
   * \f$\langle \mathbf{f}(\mathbf{r}),\,\nabla \times \overline{\overline G},\,\mathbf{f}(\mathbf{r}^{\prime})\rangle\f$.
   * @param f Testing RWG function (on observation triangle).
   * @param fp Basis RWG function (on source triangle).
   * @return Complex scalar entry of \f$\mathbf K\f$.
   * @note Not used here.
   */
  virtual dcmplx IntegrateK(RWGFun* f, RWGFun* fp);
  
  /**
   * @brief @f$\boldsymbol \delta@f$ term at an observation point.
   * @param r Observation point \f$\mathbf r\f$.
   * @param fp Source RWG function.
   * @return Complex 3-vector contribution \f$\boldsymbol \delta(\mathbf r)\f$.
   */
  virtual cvec IntegrateDelta(rvec r, RWGFun* fp);
  
  /**
   * @brief @f$\boldsymbol \kappa@f$ term at an observation point.
   * @param r Observation point \f$\mathbf r\f$.
   * @param fp Source RWG function.
   * @return Complex 3-vector contribution \f$\boldsymbol \kappa(\mathbf r)\f$.
   */
  virtual cvec IntegrateKappa(rvec r, RWGFun* fp);

  /**
   * @brief Pair of D elements with swapped triangle roles.
   * @param f Testing RWG function.
   * @param fp Basis RWG function.
   * @return \f$( D(f,fp),\, D(fp,f) )\f$.
   * @note Returns the same value for GreenFLayered3D.
   */
  virtual std::pair<dcmplx, dcmplx> PairD(RWGFun* f, RWGFun* fp);
  
  /**
   * @brief Pair of K elements with swapped triangle roles.
   * @param f Testing RWG function.
   * @param fp Basis RWG function.
   * @return \f$( K(f,fp),\, K(fp,f) )\f$.
   * @note Returns the same value for GreenFLayered3D.
   */
  virtual std::pair<dcmplx, dcmplx> PairK(RWGFun* f, RWGFun* fp);

  // -------------------------- Batched assembly helpers --------------------------
  /**
   * @brief Triangle–triangle D/K blocks assembly for all RWG pairs on two triangles.
   * @param fvec RWGs on test triangle T.
   * @param fpvec RWGs on basis triangle Tp.
   * @param DK Output matrix with pairs {D,K} for every (f,fp) pair.
   * @details Chooses between free-space and layered contributions; for layered
   * media, splits smooth (tabulated) and quasistatic parts and adds
   * optional interface line-integrals when needed.
   */
  virtual void SameTriDK(const std::vector<RWGFun*>& fvec, 
                         const std::vector<RWGFun*>& fpvec,
                         std::vector<std::vector<std::pair<dcmplx, dcmplx>>>& DK);
  /**
   * @brief Same as @ref SameTriDK and additionally accumulates gradients.
   * @param fvec RWGs on the observation triangle.
   * @param fpvec RWGs on the source triangle.
   * @param DK Output matrix of \f$(D,K)\f$ pairs.
   * @param DKG Output matrix of gradient pairs per entry.
   * @note Not implemented for layered media.
   */
  virtual void SameTriDKG(const std::vector<RWGFun*>& fvec, 
                          const std::vector<RWGFun*>& fpvec,
                          std::vector<std::vector<std::pair<dcmplx, dcmplx>>>& DK,
                          std::vector<std::vector<std::pair<cvec, cvec>>>& DKG);

  /**
   * @brief Point–triangle delta/kappa contributions for RHS assembly.
   * @param r Observation point.
   * @param fpvec RWG list on the source triangle.
   * @param DeltaKappa Output pair per RWG: (delta vector, kappa vector).
   * @details Uses near-singular smoothing and, for layered media, augments with
   * reflected/quasistatic parts from the tables.
   */
  virtual void SameTriDeltaKappa(rvec r, const std::vector<RWGFun*>& fpvec,
                                 std::vector<std::pair<cvec, cvec>>& DeltaKappa);

  // ----------------------------------------------------------------------------------                                 
  // -------------------------- Tabulation and interpolation --------------------------
  // ----------------------------------------------------------------------------------
  /**
   * @brief Fresnel reflection/transmission at interface (delegates to utils).
   * @param layerIndexI First layer index.
   * @param layerIndexJ Second layer index.
   * @param krho Transverse wavenumber.
   * @param waves "TE" or "TM".
   * @return \f$\{r, t\}\f$ complex coefficients.
   */
  std::vector<dcmplx> FresnelCoeff(int layerIndexI, int layerIndexJ, 
                                   dcmplx krho, const std::string& waves) const;
  
  /**
   * @brief Fill intra-layer tables for bottom/top half-spaces.
   * @param rq_te Quasi-static TE reflection coefficient at \f$k_{\rho}\f$ -> 0.
   * @param rq_tm Quasi-static TM reflection coefficient at \f$k_{\rho}\f$ -> 0.
   * @param inGrid Interpolation grid for this bucket.
   */
  void FillIntraTopBottom(const dcmplx& rq_te, const dcmplx& rq_tm, const Grid& inGrid);
  
  /**
   * @brief Fill intra-layer tables for internal layers (plus/minus sides).
   * @param rq_te \f$\{r_q^- , r_q^+\}\f$ TE quasi-static coefficients at \f$k_{\rho}\f$ -> 0.
   * @param rq_tm \f$\{r_q^- , r_q^+\}\f$ TM quasi-static coefficients at \f$k_{\rho}\f$ -> 0.
   * @param inGrid Interpolation grid for this bucket.
   */
  void FillIntraInside(const std::vector<dcmplx>& rq_te, 
                       const std::vector<dcmplx>& rq_tm, const Grid& inGrid);
        
  /**
   * @brief Fill inter-layer tables between two different layers.
   * @param tq_te Quasi-static TE transmission coefficient at \f$k_{\rho}\f$ -> 0.
   * @param tq_tm Quasi-static TM transmission coefficient at \f$k_{\rho}\f$ -> 0.
   * @param inGrid 3D grid \f$(r, z1, z2)\f$ for this pair.
   */
  void FillInter(const dcmplx& tq_te, const dcmplx& tq_tm, const Grid& inGrid);

  /**
   * @brief Fill all interpolation tables with Green's function values.
   * @details Loops over all interaction buckets and calls the respective
   * filling routines.
   */
  void FillGreenFTable();

  /**
   * @brief Test filling of interpolation tables with Green's function values.
   * @param intraTBType Intra-layer top/bottom selector.
   * @param intraPMType Intra-layer plus/minus selector.
   * @param interType Inter-layer coefficient family selector.
   * @param numPoints Number of sample points to generate per table/dimension.
   */
  void testFillGreenFTable(const std::string& intraTBType, const std::string& intraPMType, 
                           const std::string& interType, const int& numPoints);
                        
   /**
   * @brief Evaluate 2D interpolant at one query point.
   * @param x Query X (e.g., radial distance r).
   * @param y Query Y (e.g., z-like coordinate).
   * @param inCoeffs Tuple produced by tabulation: (F_flat, px, py, Xc, Yc).
   * F_flat is column-major (x fastest) complex values on grid.
   * @param type Which component of F to interpolate: "real" or "imag".
   * @return Interpolated value (real-valued) of the requested component.
   * @note Uses regridpack (bicubic by default). Inputs are clamped to grid box.
   * @warning If `ier!=0` is reported, the regridpack call failed; result is still
   * returned but may be inaccurate.
   */
  double evaluateSinglePoint(double x, double y, 
                             std::tuple<std::vector<dcmplx>, std::vector<int>, std::vector<int>, 
                                        std::vector<double>, std::vector<double>>& inCoeffs, 
                             std::string type);

  /**
   * @brief Evaluate 3D interpolant at one query point using per-slice 2D regridding
   * and Catmull–Rom cubic in z.
   * @param x Query X (typically r).
   * @param y Query Y (typically z1).
   * @param z Query Z (typically z2).
   * @param inCoeffs Tuple (Xv,Yv,Zv,Fc) where Fc is complex column-major volume.
   * @param type "real" or "imag" to select part of Fc.
   * @return Interpolated value (real-valued) at (x,y,z).
   * @details Performs bicubic 2D interpolation on bracketing z-slices and blends 
   * along z with Catmull–Rom. Clamps queries to the grid box; uses linear interpolation 
   * if mz==1.
   */
  double evaluateSinglePoint(
                        double x, double y, double z, 
                        std::tuple<std::vector<double>, std::vector<double>, std::vector<double>, 
                                    std::vector<dcmplx>>& inCoeffs, 
                        std::string type);

  /**
   * @brief Write complex data with associated 2D coordinates to a text file.
   * @param x X-coordinates (size N).
   * @param y Y-coordinates (size N).
   * @param inVec Complex values at (x,y) (size N).
   * @param filename Output .txt filename (overwritten if exists).
   */
  void writeVectorToFile(const std::vector<double>& x, const std::vector<double>& y, 
                         const std::vector<dcmplx>& inVec, const std::string& filename);
  
  /**
   * @brief Write real data with associated 2D coordinates (rho as rvec) to a text file.
   * @param rho 2D positions as rvec entries (size N); only in-plane components are used.
   * @param inVec Real values at positions (size N).
   * @param filename Output .txt filename (overwritten if exists).
   */                         
  void writeVectorToFile(const std::vector<rvec>& rho, const std::vector<double>& inVec, 
                         const std::string& filename);
  
  /**
   * @brief Write complex data with associated 3D coordinates to a text file.
   * @param x X-coordinates (size N).
   * @param y Y-coordinates (size N).
   * @param z Z-coordinates (size N).
   * @param inVec Complex values at (x,y,z) (size N).
   * @param filename Output .txt filename (overwritten if exists).
   */
  void writeVectorToFile(const std::vector<double>& x, const std::vector<double>& y, 
                         const std::vector<double>& z, const std::vector<dcmplx>& inVec, 
                         const std::string& filename);

  /**
   * @brief Evaluate Green's function components for **intra-layer** interactions by
   * interpolating the smooth (Sommerfeld) coefficient tables. Interpolates pretabulated 
   * complex coefficients on the (R, Z1, Z2) grid built during @p FillGreenFTable for 
   * the requested bucket.
   * @param interactionType Identifier for the intra-layer type. Must match the key used 
   * at tabulation time.
   * @param R In-plane radial separation @f$\rho@f$ between source and observation.
   * @param Z1 Vertical distance z+z' for half-space intra-layer interactions or z-z'
   * for inner intra-layer interactions (the different definitions, can be found 
   * in @p cartesianToLayeredTab() of @p LayeredMediaUtils).
   * @param Z2 Vertical distance z+z' inner intra-layer (the definition can be found 
   * in @p cartesianToLayeredTab() of @p LayeredMediaUtils).
   * @param coeffIndex Index selecting which interpolation coefficients to use.
   * @return std::map from component name to complex value at (R, Z1, Z2). The key set 
   * follows the tabulation.
   * @note Queries are clamped to the grid box; outside-range inputs fall back to
   * boundary interpolation.
   */
  std::map<std::string,dcmplx> evalGreenSmoothIntra(const std::string& interactionType,
                                                    double R, double Z1, double Z2, 
                                                    int coeffIndex);

  /**
   * @brief Evaluate Green's-function components for **inter-layer** interactions by
   * interpolating the smooth (Sommerfeld) coefficient volume. Uses the 3D tables built 
   * for pairs of distinct layers; indices and table keys must match those produced in 
   * @p FillGreenFTable.
   * @param R In-plane radial separation @f$\rho@f$ between source and observation.
   * @param Z1 Vertical distance z inter-layer interaction (definition can 
   * be found in @p cartesianToLayeredTab() of @p LayeredMediaUtils).
   * @param Z2 Vertical distance z' inter-layer (definition can be found 
   * in @p cartesianToLayeredTab() of @p LayeredMediaUtils).
   * @param coeffIndex Index selecting which interpolation coefficients to use.
   * @return std::map from component name to complex value at (R, Z1, Z2). The key set 
   * follows the tabulation.
   * @note Queries are clamped to the volume bounds; linear blending may be used
   * along sparse axes depending on grid density.
   */
  std::map<std::string,dcmplx> evalGreenSmoothInter(double R, double Z1, double Z2, 
                                                    int coeffIndex);

  /**
   * @brief Evaluate **quasi-static (QS)** contributions for **intra-layer** interactions.
   * @param pos Reference Cartesian position used to determine the type of intra-layer
   * interaction and the appropriate QS images.
   * @param interactionType Interaction identifier (e.g., "intraTB" or "intraIn").
   * @param rho In-plane separation @f$\rho@f$ between the two points.
   * @param Z1 Vertical distance z+z' for half-space intra-layer interactions or z-z'
   * for inner intra-layer interactions (the different definitions, can be found 
   * in @p cartesianToLayeredQS() of @p LayeredMediaUtils).
   * @param Z2 Vertical distance z+z' inner intra-layer (the definition can be found 
   * in @p cartesianToLayeredQS() of @p LayeredMediaUtils).
   * @param R1 Full 3D distance in (rho, Z1).
   * @param R2 Full 3D distance in (rho, Z2).
   * @param singCancFields If true, arrange/return terms for singular-cancellation
   * formulations used by line/surface splits.
   * @return Map from component name to complex QS value. The key set follows the tabulation.
   */
  std::map<std::string, dcmplx> evalGreenQSIntra(const rvec pos, const std::string& interactionType,
                                                 double rho, double Z1, double Z2, double R1, double R2,
                                                 bool singCancFields = false);
                                                
  /**
   * @brief Evaluate **quasi-static (QS)** contributions for **inter-layer** interactions.
   * @param pos1 Observation point (Cartesian) in layer i.
   * @param pos2 Source point (Cartesian) in layer j (j != i).
   * @param rho In-plane separation @f$\rho@f$ between pos1 and pos2.
   * @param Z Vertical distance between the two points.
   * @param R Full 3D source–observer separation.
   * @param singCancFields If true, arrange/return terms for singular-cancellation
   * formulations used by line/surface splits.
   * @return Map from component name to complex QS value. The key set follows the tabulation.
   * @note Uses the low-krho limits of the Fresnel coefficients as implemented in
   * the layered-media utilities.
   */
  std::map<std::string, dcmplx> evalGreenQSInter(const rvec pos1, const rvec pos2, 
                                                 double rho, double Z, double R,
                                                 bool singCancFields = false);

  /**
   * @brief Test assembly of same-triangle D and K matrix elements.
   * @param verts1 Vertices of the first triangle (3 entries: r0,r1,r2).
   * @param verts2 Vertices of the second triangle (3 entries: r0',r1',r2').
   * @param gType Green's function tag used in the test.
   * @param vertPrint Verbose label or mode controlling vertex printing in the debug output.
   */
  void testSameTriDK(std::vector<rvec> verts1, std::vector<rvec> verts2, 
                     std::string gType, std::string vertPrint);
};

#endif
