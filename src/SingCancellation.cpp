/**
 * @file SingCancellation.cpp
 * @brief Implementation of the singularity-cancellation paired quadrature builders.
 */

#include "SingCancellation.h"

// Construct and immediately build the paired quadrature.
SingCancellation::SingCancellation(int commonNodes, Triangle* tPtr1, Triangle* tPtr2,
                                   int nQuadrature, double(*xQuadrature)[2], double* wQuadrature)
    : commonNodes(commonNodes), tPtr1(tPtr1), tPtr2(tPtr2), nQuadrature(nQuadrature),
      xQuadrature(xQuadrature), wQuadrature(wQuadrature) {
        if (commonNodes == 1) commonVertex();
        else if (commonNodes == 2) commonEdge();
        else if (commonNodes == 3) commonFacet();
        else {
          throw std::invalid_argument("commonNodes must be 1, 2, or 3.");
        }
      }

// Trivial destructor.
SingCancellation::~SingCancellation() { }

/**
 * @details Build quadrature for triangles with one *common vertex*.
 * Creates the raw 4D tensor-product rule (@f$ n^4 @f$ points), then applies
 * the two-way symmetrization (@f$ \xi \leftrightarrow \eta @f$) to obtain 
 * @f$ 2 n^4 @f$ paired samples. Each sample produces one barycentric triple 
 * on triangle #1 and the corresponding triple on triangle #2, plus the combined weight:
 * @f$ w = 4 w_i w_j w_p w_q J(\xi, \eta, ...) @f$.
 *
 * Internally uses Triangle::FindCommonVertexIndices() to determine how the
 * barycentric coordinates must be permuted on each triangle so that the
 * *same* physical vertex is aligned in both local parameterizations.
 *
 * @warning On failure to detect the shared vertex, a runtime_error is thrown.
 */
void SingCancellation::commonVertex() {
  int size = nQuadrature * nQuadrature * nQuadrature * nQuadrature; // number of raw 4D points

  // allocate working arrays:
  //    xi1, xi2     : size x 2 (1 + 1 symmetric)
  //    eta1,eta2    : size x 2
  //    jac, weights : size
  std::vector<std::vector<double>> xi1(size, std::vector<double>(2));
  std::vector<std::vector<double>> xi2(size, std::vector<double>(2));
  std::vector<std::vector<double>> eta1(size, std::vector<double>(2));
  std::vector<std::vector<double>> eta2(size, std::vector<double>(2));
  std::vector<double> jac(size), weights(size);

  // Build the raw tensor-product nodes/weights
  for(int i = 0, idx = 0; i < nQuadrature; ++i) {
    for(int j = 0; j < nQuadrature; ++j) {
      for(int p = 0; p < nQuadrature; ++p) {
        for(int q = 0; q < nQuadrature; ++q, ++idx) {
          xi1[idx][0]  = xQuadrature[q][0];
          xi2[idx][0]  = xQuadrature[q][0] * xQuadrature[p][0];
          eta1[idx][0] = xQuadrature[q][0] * xQuadrature[j][0];
          eta2[idx][0] = xQuadrature[q][0] * xQuadrature[j][0] * xQuadrature[i][0];

          // Combined Jacobian for the nested maps
          jac[idx]  = xQuadrature[j][0] * xQuadrature[q][0] * xQuadrature[q][0] * xQuadrature[q][0];
          
          // 4D tensor-product weight
          weights[idx] = wQuadrature[q] * wQuadrature[p] * wQuadrature[j] * wQuadrature[i];
        }
      }
    }
  }
  
  for(int i = 0; i < size; ++i)
  {
    // --- symmetrize: append (ξ->η, η->ξ) to get 12 total
    // 1st half of sym: η <- ξ
    eta1[i][1] = xi1[i][0];
    eta2[i][1] = xi2[i][0];
    // 2nd half: ξ <- η
    xi1[i][1] = eta1[i][0];
    xi2[i][1] = eta2[i][0];
  }

  // Output arrays
  pointsMap1.assign( 2 * size, std::vector<double>( 3, 0.0 ) );
  pointsMap2.assign( 2 * size, std::vector<double>( 3, 0.0 ) );
  w.resize(2 * size);

  // Determine alignment/permutation on each triangle
  std::pair<int,int> commonVertex = tPtr1->FindCommonVertexIndices(tPtr2);
  if (commonVertex.first == 0 && commonVertex.second == 0) {
    throw std::runtime_error("No common vertex found between the triangles.");
  }
  std::vector<int> perm1 = this->getPermutation(commonVertex.first);
  std::vector<int> perm2 = this->getPermutation(commonVertex.second);

  // Fill paired barycentric nodes and weights
  for(int i = 0; i < size; ++i) {
    for(int k = 0; k < 2; ++k) {
      int j = k * size + i; // index in the 1D vector
      pointsMap1[j][perm1[0]] = 1.0 - xi1[i][k];
      pointsMap1[j][perm1[1]] = xi1[i][k] - xi2[i][k];
      pointsMap1[j][perm1[2]] = xi2[i][k];
      pointsMap2[j][perm2[0]] = 1.0 - eta1[i][k];
      pointsMap2[j][perm2[1]] = eta1[i][k] - eta2[i][k];
      pointsMap2[j][perm2[2]] = eta2[i][k];
      w[j]  = 4.0 * weights[i] * jac[i];
    }
  }
#ifdef DEBUG
    // Write out five MATLAB-formatted .txt files:
    writeMatlabMatrixToTxt("x1.txt", "x1_cpp", "x", pointsMap1, 2 * size, 6);
    writeMatlabMatrixToTxt("y1.txt", "y1_cpp", "y", pointsMap1, 2 * size, 6);
    writeMatlabMatrixToTxt("z1.txt", "z1_cpp", "z", pointsMap1, 2 * size, 6);
    writeMatlabMatrixToTxt("x2.txt", "x2_cpp", "x", pointsMap2, 2 * size, 6);
    writeMatlabMatrixToTxt("y2.txt", "y2_cpp", "y", pointsMap2, 2 * size, 6);
    writeMatlabMatrixToTxt("z2.txt", "z2_cpp", "z", pointsMap2, 2 * size, 6);
    writeMatlabVectorToTxt("w.txt",  "w_cpp",  w, 6);
#endif
}

/**
 * @details Build quadrature for triangles sharing a *common edge*.
 * Mapping creates six base cases (@f$k=0...5@f$), then symmetrizes to twelve,
 * yielding @f$12n^4@f$ paired nodes overall. The @f$\omega@f$-dependent Jacobian
 * @f$ x_1 \omega^2 (1-\omega) @f$ is folded into the weights.
 *
 * @warning std::runtime_error if no common edge is found.
 */
void SingCancellation::commonEdge() {
  int size = nQuadrature * nQuadrature * nQuadrature * nQuadrature; // number of raw 4D points

  // build the raw 4-D tensor product grid of (ω, x1, x2, χ1) and weights
  std::vector<double> omega(size), x1(size), x2(size), chi1(size), weights(size);
  for(int i = 0, idx = 0; i < nQuadrature; ++i) {
    for(int j = 0; j < nQuadrature; ++j) {
      for(int p = 0; p < nQuadrature; ++p) {
        for(int q = 0; q < nQuadrature; ++q, ++idx) {
          omega[idx] = xQuadrature[q][0];
          x1[idx] = xQuadrature[p][0];
          x2[idx] = xQuadrature[j][0];
          chi1[idx] = xQuadrature[i][0];
          weights[idx] = wQuadrature[q] * wQuadrature[p] * wQuadrature[j] * wQuadrature[i];
        }
      }
    }
  }

  // allocate working arrays:
  //    mu1, mu2  : size x 6
  //    xi1, xi2  : size x 12 (6 + 6 symmetric)
  //    eta1,eta2 : size x 12
  //    jac       : size
  std::vector<std::vector<double>> mu1(size, std::vector<double>(6));
  std::vector<std::vector<double>> mu2(size, std::vector<double>(6));
  std::vector<std::vector<double>> xi1(size, std::vector<double>(12));
  std::vector<std::vector<double>> xi2(size, std::vector<double>(12));
  std::vector<std::vector<double>> eta1(size, std::vector<double>(12));
  std::vector<std::vector<double>> eta2(size, std::vector<double>(12));
  std::vector<double> jac(size);
  for (int i = 0; i < size; ++i) {
    // 6 base cases (k = 0..5)
    mu1[i][0] = - omega[i] * x1[i];
    mu2[i][0] = - omega[i] * x1[i] * x2[i];
    xi1[i][0] = ( 1.0 - omega[i] ) * chi1[i] + omega[i];
    xi2[i][0] = omega[i] * ( 1.0 - x1[i] + x1[i] * x2[i] );
    
    mu1[i][1] = omega[i] * x1[i];
    mu2[i][1] = omega[i] * x1[i] * x2[i];
    xi1[i][1] = ( 1.0 - omega[i] ) * chi1[i] + omega[i] * ( 1.0 - x1[i] );
    xi2[i][1] = omega[i] * ( 1.0 - x1[i] );
    
    mu1[i][2] = - omega[i] * x1[i] * x2[i];
    mu2[i][2] = omega[i] * x1[i] * ( 1.0 - x2[i] );
    xi1[i][2] = ( 1.0 - omega[i] ) * chi1[i] + omega[i];
    xi2[i][2] = omega[i] * ( 1.0 - x1[i] );
    
    mu1[i][3] = omega[i] * x1[i] * x2[i];
    mu2[i][3] = omega[i] * x1[i] * ( x2[i] - 1.0 );
    xi1[i][3] = ( 1.0 - omega[i] ) * chi1[i] + omega[i] * ( 1.0 - x1[i] * x2[i] );
    xi2[i][3] = omega[i] * ( 1.0 - x1[i] * x2[i] );
    
    mu1[i][4] = - omega[i] * x1[i] * x2[i];
    mu2[i][4] = - omega[i] * x1[i];
    xi1[i][4] = ( 1.0 - omega[i] ) * chi1[i] + omega[i];
    xi2[i][4] = omega[i];
    
    mu1[i][5] = omega[i] * x1[i] * x2[i];
    mu2[i][5] = omega[i] * x1[i];
    xi1[i][5] = ( 1.0 - omega[i] ) * chi1[i] + omega[i] * ( 1.0 - x1[i] * x2[i] );
    xi2[i][5] = omega[i] * ( 1.0 - x1[i] );

    jac[i] = x1[i] * omega[i] * omega[i] * ( 1.0 - omega[i] );

    for(int k = 0; k < 6; ++k)
    {
      // --- compute η = μ + ξ for the 6 cases
      eta1[i][k] = mu1[i][k] + xi1[i][k];
      eta2[i][k] = mu2[i][k] + xi2[i][k];

      // --- symmetrize: append (ξ->η, η->ξ) to get 12 total
      // 1st half of sym: η <- ξ
      eta1[i][k + 6] = xi1[i][k];
      eta2[i][k + 6] = xi2[i][k];
      // 2nd half: ξ <- η
      xi1[i][k + 6] = eta1[i][k];
      xi2[i][k + 6] = eta2[i][k];
    }
  }

  pointsMap1.assign( 12 * size, std::vector<double>( 3, 0.0 ) );
  pointsMap2.assign( 12 * size, std::vector<double>( 3, 0.0 ) );
  w.resize(12 * size);

  // Determine the shared edge and permutations
  std::pair<int,int> commonEdges = tPtr1->FindCommonEdgeIndices(tPtr2);
  if (commonEdges.first == 0 && commonEdges.second == 0) {
    throw std::runtime_error("No common edge found between the triangles.");
  }
  std::vector<int> perm1 = this->getPermutation(commonEdges.first);
  std::vector<int> perm2 = this->getPermutation(commonEdges.second);
  
  // Fill the points and weights for both maps
  for(int i = 0; i < size; ++i) {
    for(int k = 0; k < 12; ++k) {
      int j = k * size + i; // index in the 1D vector
      pointsMap1[j][perm1[0]] = 1.0 - xi1[i][k];
      pointsMap1[j][perm1[1]] = xi1[i][k] - xi2[i][k];
      pointsMap1[j][perm1[2]] = xi2[i][k];
      pointsMap2[j][perm2[0]] = 1.0 - eta1[i][k];
      pointsMap2[j][perm2[1]] = eta1[i][k] - eta2[i][k];
      pointsMap2[j][perm2[2]] = eta2[i][k];
      w[j]  = 0.5 * (4.0 * weights[i] * jac[i]);  // factor 1/2 from mapping
    }
  }
#ifdef DEBUG
    // Write out five MATLAB-formatted .txt files:
    writeMatlabMatrixToTxt("x1.txt", "x1_cpp", "x", pointsMap1, 12 * size, 6);
    writeMatlabMatrixToTxt("y1.txt", "y1_cpp", "y", pointsMap1, 12 * size, 6);
    writeMatlabMatrixToTxt("z1.txt", "z1_cpp", "z", pointsMap1, 12 * size, 6);
    writeMatlabMatrixToTxt("x2.txt", "x2_cpp", "x", pointsMap2, 12 * size, 6);
    writeMatlabMatrixToTxt("y2.txt", "y2_cpp", "y", pointsMap2, 12 * size, 6);
    writeMatlabMatrixToTxt("z2.txt", "z2_cpp", "z", pointsMap2, 12 * size, 6);
    writeMatlabVectorToTxt("w.txt",  "w_cpp",  w, 6);
#endif
}

/**
 * @details Build quadrature for *identical* triangles (common facet).
 * Uses three base mappings and symmetry (@f$ \xi \leftrightarrow \eta @f$) 
 * to produce @f$ 6 n^4 @f$ pairs. A per-case Jacobian (three values) is 
 * applied to the weights.
 */
void SingCancellation::commonFacet() {
  int size = nQuadrature * nQuadrature * nQuadrature * nQuadrature; // number of raw 4D points

  // build the raw 4-D tensor product grid of (ω, x1, x2, χ1) and weights
  std::vector<double> omega(size), x(size), chi1(size), chi2(size), weights(size);
  for(int i = 0, idx = 0; i < nQuadrature; ++i) {
    for(int j = 0; j < nQuadrature; ++j) {
      for(int p = 0; p < nQuadrature; ++p) {
        for(int q = 0; q < nQuadrature; ++q, ++idx) {
          omega[idx] = xQuadrature[q][0];
          x[idx] = xQuadrature[p][0];
          chi1[idx] = xQuadrature[j][0];
          chi2[idx] = xQuadrature[i][0];
          weights[idx] = wQuadrature[q] * wQuadrature[p] * wQuadrature[j] * wQuadrature[i];
        }
      }
    }
  }

  // allocate working arrays:
  //    mu1, mu2  : size x 3
  //    xi1, xi2  : size x 6 (3 + 3 symmetric)
  //    eta1,eta2 : size x 6
  //    jac       : size x 6
  std::vector<std::vector<double>> mu1(size, std::vector<double>(3));
  std::vector<std::vector<double>> mu2(size, std::vector<double>(3));
  std::vector<std::vector<double>> xi1(size, std::vector<double>(6));
  std::vector<std::vector<double>> xi2(size, std::vector<double>(6));
  std::vector<std::vector<double>> eta1(size, std::vector<double>(6));
  std::vector<std::vector<double>> eta2(size, std::vector<double>(6));
  std::vector<std::vector<double>> jac(size, std::vector<double>(6));

  for (int i = 0; i < size; ++i) {
    // k = 0
    mu1[i][0] = omega[i];
    mu2[i][0] = omega[i] * x[i];
    xi1[i][0] = ( 1.0 - mu1[i][0] ) * chi1[i];
    xi2[i][0] = xi1[i][0] * chi2[i];
    jac[i][0] = omega[i] * ( 1.0 - mu1[i][0] ) * xi1[i][0];
    
    // k = 1
    mu1[i][1] = omega[i] * x[i];
    mu2[i][1] = omega[i] * ( x[i] - 1.0 );
    xi1[i][1] = ( 1.0 - mu1[i][1] + mu2[i][1] ) * chi1[i] - mu2[i][1];
    xi2[i][1] = ( xi1[i][1] + mu2[i][1] ) * chi2[i] - mu2[i][1];
    jac[i][1] = omega[i] * ( 1.0 - mu1[i][1] + mu2[i][1] ) * ( xi1[i][1] + mu2[i][1] );
    
    // k = 2
    mu1[i][2] = omega[i] * x[i];
    mu2[i][2] = omega[i];
    xi1[i][2] = ( 1.0 - mu2[i][2] ) * chi1[i] + mu2[i][2] - mu1[i][2];
    xi2[i][2] = ( xi1[i][2] - mu2[i][2] + mu1[i][2] ) * chi2[i];
    jac[i][2] = omega[i] * ( 1.0 - mu2[i][2] ) * ( xi1[i][2] - mu2[i][2] + mu1[i][2] );

    for(int k = 0; k < 3; ++k)
    {
      // --- compute η = μ + ξ for the 3 cases
      eta1[i][k] = mu1[i][k] + xi1[i][k];
      eta2[i][k] = mu2[i][k] + xi2[i][k];

      // --- symmetrize: append (ξ->η, η->ξ) to get 12 total
      // 1st half of sym: η <- ξ
      eta1[i][k + 3] = xi1[i][k];
      eta2[i][k + 3] = xi2[i][k];
      // 2nd half: ξ <- η
      xi1[i][k + 3] = eta1[i][k];
      xi2[i][k + 3] = eta2[i][k];
    }
  }

  pointsMap1.assign( 6 * size, std::vector<double>( 3, 0.0 ) );
  pointsMap2.assign( 6 * size, std::vector<double>( 3, 0.0 ) );
  w.resize(6 * size);

  // Fill paired nodes and weights
  for(int i = 0; i < size; ++i) {
    for(int k = 0; k < 6; ++k) {
      int j = k * size + i; // index in the 1D vector
      pointsMap1[j][0] = 1.0 - xi1[i][k];
      pointsMap1[j][1] = xi1[i][k] - xi2[i][k];
      pointsMap1[j][2] = xi2[i][k];
      pointsMap2[j][0] = 1.0 - eta1[i][k];
      pointsMap2[j][1] = eta1[i][k] - eta2[i][k];
      pointsMap2[j][2] = eta2[i][k];
      w[j]  = 4.0 * weights[i] * jac[i][k % 3];
    }
  }
#ifdef DEBUG
    // Write out five MATLAB-formatted .txt files:
    writeMatlabMatrixToTxt("x1.txt", "x1_cpp", "x", pointsMap1, 6 * size, 6);
    writeMatlabMatrixToTxt("y1.txt", "y1_cpp", "y", pointsMap1, 6 * size, 6);
    writeMatlabMatrixToTxt("z1.txt", "z1_cpp", "z", pointsMap1, 6 * size, 6);
    writeMatlabMatrixToTxt("x2.txt", "x2_cpp", "x", pointsMap2, 6 * size, 6);
    writeMatlabMatrixToTxt("y2.txt", "y2_cpp", "y", pointsMap2, 6 * size, 6);
    writeMatlabMatrixToTxt("z2.txt", "z2_cpp", "z", pointsMap2, 6 * size, 6);
    writeMatlabVectorToTxt("w.txt",  "w_cpp",  w, 6);
#endif
}

/**
 * @details Map "shift code" -> permutation indices of a barycentric triple.
 * Encodes the rotations/reflections needed to align local 
 * @f$ (\lambda_1, \lambda_2, \lambda_3) @f$ choices across the two triangles 
 * so the *same physical* vertex/edge is addressed.
 */
std::vector<int> SingCancellation::getPermutation(int shift) {
  switch (shift) {
    case  1: // [x,y] = (x,y)
      return {0, 1, 2};
    case  2: // [x,y] = (z,x)
      return {2, 0, 1};
    case  3: // [x,y] = (y,z)
      return {1, 2, 0};
    case -1: // [x,y] = (x,z)
      return {0, 2, 1};
    case -2: // [x,y] = (z,y)
      return {2, 1, 0};
    case -3: // [x,y] = (y,x)
      return {1, 0, 2};
    default:
      throw std::invalid_argument("getPermutation: unexpected shift");
  }
}

//Dump one barycentric component (x/y/z) as a MATLAB column vector.
void SingCancellation::writeMatlabMatrixToTxt(
    const std::string &filename, const std::string &varname, const std::string &component,
    const std::vector<std::vector<double>> &V, int size, int prec)
{
  std::ofstream fout(filename);
  if(!fout) {
      std::cerr << "Error: could not open " 
                << filename << " for writing\n";
      return;
  }

  int comp = -1; // Default component index
  if (component == "x") {
    comp = 0; // x component
  } else if (component == "y") {
    comp = 1; // y component
  } else if (component == "z") {
    comp = 2; // z component
  } else if (component != "") {
    std::cerr << "Error: Invalid component '" << component << "' specified.\n";
    return;
  }

  // Start the MATLAB assignment
  fout << varname << " = [\n"
        << std::fixed << std::setprecision(prec);

  // Write each row
  for(int i = 0; i < size; ++i) {
    fout << V[i][comp] << "; ";
  }
  fout << "];\n";
  fout.close();
}

//Dump a real vector as a MATLAB column vector.
void SingCancellation::writeMatlabVectorToTxt(
    const std::string &filename,
    const std::string &varname,
    const std::vector<double> &w,
    int prec)
{
    std::ofstream fout(filename);
    if (!fout) {
        std::cerr << "Error: could not open "
                  << filename << " for writing\n";
        return;
    }

    // MATLAB assignment
    fout << varname << " = [\n"
         << std::fixed << std::setprecision(prec);

    for (double v : w) {
        fout << v << ";\n";
    }

    fout << "];\n";
}
