/**
 * @file SingCancellation.h
 * @brief Quadrature point/weight builders for singularity-cancellation on triangles.
 *
 * This utility constructs *paired* quadrature rules on two triangles that
 * share a **vertex**, an **edge**, or have the **same facet**. The rules follow
 * mappings used in Taylor (2003), D'Elia et al. (2011), Sarraf et al. (2014),
 * and produce synchronized barycentric nodes on each triangle together with
 * compensated Jacobians to cancel kernel singularities.
 *
 * Usage pattern:
 *   - Construct with the pair of triangles and a chosen 1D Gauss–Legendre rule
 *      (n = 1, 4, or 8 nodes is typical: see quadrature.{h,cpp}).
 *   - Internally, the constructor detects how the triangles touch (1,2,3 common nodes)
 *      and calls the matching builder: commonVertex(), commonEdge(), or commonFacet().
 *   - Retrieve the paired barycentric points via getPointsMap() and the weights via getWeights().
 */

#ifndef SINGCANCELLATION_H
#define SINGCANCELLATION_H

#include "Triangle.h"
#include "RWGFun.h"
#include "quadrature.h"

/// \ingroup quadrature
/**
 * @class SingCancellation
 * @brief Builds synchronized quadrature nodes/weights for a pair of touching triangles.
 *
 * The class generates *paired* sampling points (barycentric) on the two triangles
 * together with properly scaled weights so that singular and near-singular integrals
 * (e.g., double-surface integrals with 1/R kernels) are integrable with high accuracy.
 *
 * After construction, call:
 *   - getPointsMap(0) for triangle #1 nodes
 *   - getPointsMap(1) for triangle #2 nodes
 *   - getWeights() for the weights associated with each *paired* sample
 */
class SingCancellation
{
private:
    int commonNodes;            ///< Number of common nodes between two triangles
    Triangle* tPtr1;            ///< First triangle pointer
    Triangle* tPtr2;            ///< Second triangle pointer
    int nQuadrature;            ///< Number of quadrature points for integration
    double(*xQuadrature)[2];    ///< Quadrature points for integration
    double* wQuadrature;        ///< Weights for quadrature points
    std::vector<double> w;      ///< Weights for integration
    
    /// Data structure to hold quadrature points and weights for the pair of triangles
    std::vector<std::vector<double>> pointsMap1;
    /// Data structure to hold quadrature points and weights for the pair of triangles
    std::vector<std::vector<double>> pointsMap2;
    
public:
    /**
     * @brief Construct and immediately generate the paired quadrature for a triangle pair.
     *
     * Detects how the triangles touch by @p commonNodes and calls the appropriate builder:
     *  - 1 -> commonVertex()
     *  - 2 -> commonEdge()
     *  - 3 -> commonFacet()
     *
     * @param commonNodes   Number of common nodes (1, 2, or 3).
     * @param tPtr1         Pointer to first triangle.
     * @param tPtr2         Pointer to second triangle.
     * @param nQuadrature   Number of 1D GL points (typically 1, 4, or 8).
     * @param xQuadrature   GL nodes array of shape [nQuadrature][2] with @f$ (\xi, 1-\xi) @f$.
     * @param wQuadrature   GL weights array of length nQuadrature.
     *
     * @throws std::invalid_argument if @p commonNodes not in {1,2,3}.
     */
    SingCancellation(int commonNodes, Triangle* tPtr1, Triangle* tPtr2,
                     int nQuadrature, double(*xQuadrature)[2], double* wQuadrature);
    
    /** @brief Trivial destructor. */
    ~SingCancellation();
    
    /**
     * @brief Build paired quadrature for triangles sharing a *single vertex*.
     *
     * Produces @f$ 2 n^4 @f$ paired nodes (symmetrization doubles the base 
     * @f$ n^4 @f$ set), with weights including the appropriate Jacobian factors 
     * arising from the four nested scalar transforms. Based on Taylor (2003), 
     * D'Elia et al. (2011), Sarraf et al. (2014).
     *
     * @note Called automatically by the constructor when \p commonNodes==1.
     */
    void commonVertex();
    
    /**
     * @brief Build paired quadrature for triangles sharing a *common edge*.
     *
     * Produces @f$ 12 n^4 @f$ paired nodes per the six base mappings plus their
     * symmetric counterparts (@f$ \xi \leftrightarrow \eta @f$). Includes the 
     * @f$ \omega @f$-dependent Jacobian. Based on Taylor (2003), 
     * D'Elia et al. (2011), Sarraf et al. (2014).
     *
     * @note Called automatically by the constructor when \p commonNodes==2.
     */
    void commonEdge();

    /**
     * @brief Build paired quadrature for *identical triangles* (common facet).
     *
     * Produces @f$ 6 n^4 @f$ paired nodes from three base mappings and symmetry,
     * with per-case Jacobians folded into the weights. Based on Taylor (2003), 
     * D'Elia et al. (2011), Sarraf et al. (2014).
     *
     * @note Called automatically by the constructor when \p commonNodes==3.
     */
    void commonFacet();

    /**
     * @brief Permutation of barycentric components induced by vertex/edge alignment.
     *
     * The shift encodes how the @f$ (\lambda_1, \lambda_2, \lambda_3) @f$ triple 
     * must be permuted so that the first component corresponds to the vertex 
     * opposite the "free" entity.
     *
     * @param shift Encoded shift (−3...−1, 1...3).
     * @return A permutation vector of indices in {0,1,2}.
     *
     * @throws std::invalid_argument on unsupported @p shift.
     */
    std::vector<int> getPermutation(int shift);

    /**
     * @brief Access the paired barycentric nodes for triangle i.
     *
     * @param triangleIndex 0 for the first triangle, 1 for the second.
     * @return Points map of size (paired nodes) x 3.
     *
     * @throws std::invalid_argument if @p triangleIndex not in {0,1}.
     */
    std::vector<std::vector<double>> getPointsMap(int triangleIndex) {
        if (triangleIndex == 0) {
            return pointsMap1;
        } else if (triangleIndex == 1) {
            return pointsMap2;
        } else {
            throw std::invalid_argument("Invalid triangle index. Use 1 or 2.");
        }
    }

    /**
     * @brief Access the weights associated with each paired sampling point.
     * @return Vector of size paired nodes.
     */
    std::vector<double> getWeights() { return w; }

    /**
     * @brief Write a single matrix column (chosen component) as a MATLAB style vector.
     *
     * Writes:
     *   varname = [
     *     v1;
     *     v2;
     *     ...
     *   ];
     *
     * @param filename  Output text filename.
     * @param varname   MATLAB variable name to assign.
     * @param component Which barycentric component to export: "x","y","z" (0,1,2).
     *                  If empty, writes the first component by default.
     * @param V         Matrix to export (rows correspond to points).
     * @param size      Number of rows to write (typically V.size()).
     * @param prec      Decimal precision (default 6).
     */
    void writeMatlabMatrixToTxt(const std::string &filename,
                                const std::string &varname,
                                const std::string &component, 
                                const std::vector<std::vector<double>> &V,
                                int size, int prec = 6);

    /**
     * @brief Write a real vector as a MATLAB style column vector.
     *
     * @param filename  Output text filename.
     * @param varname   MATLAB variable name to assign.
     * @param w         Data vector to write.
     * @param prec      Decimal precision (default 6).
     */
    void writeMatlabVectorToTxt(const std::string &filename, const std::string &varname,
                                const std::vector<double> &w, int prec = 6);
};

#endif