// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <eigen3/Eigen/Dense>
#include <iostream>
#include <utility>

#include "/home/rehmemk/git/cddlib/lib-src/cdd.h"
#include "/home/rehmemk/git/volume_approximation/include/volume/volume.h"

/**
 *
 * @param dim dimension
 * @param A0 the active subspaces (columnwise)
 * @param point a point in the orthogonal plane
 */
dd_MatrixPtr createHPolytopeDescription(size_t dim, Eigen::MatrixXd A0, Eigen::VectorXd point) {
  dd_MatrixPtr hMatrix = dd_CreateMatrix(2 * dim + 2 * A0.cols(), dim + 1);
  for (size_t i = 0; i < dim; i++) {
    // cube [0,1]^dim
    hMatrix->matrix[i][0][0] = 1.0;
    hMatrix->matrix[i][i + 1][0] = -1.0;
    hMatrix->matrix[dim + i][i + 1][0] = 1.0;
  }

  // cutting orthogonal to (1D) active subspace (generalize for Multi D)
  Eigen::VectorXd prod = A0.col(0).transpose() * point;
  hMatrix->matrix[2 * dim][0][0] = prod(0);
  hMatrix->matrix[2 * dim + 1][0][0] = -prod(0);

  for (size_t j = 1; j < dim + 1; j++) {
    hMatrix->matrix[2 * dim][j][0] = -A0(j - 1, 0);
    hMatrix->matrix[2 * dim + 1][j][0] = A0(j - 1, 0);
  }
  hMatrix->numbtype = dd_Real;
  hMatrix->representation = dd_Inequality;
  return hMatrix;
}

Eigen::MatrixXd affineHull(dd_MatrixPtr vRep, Eigen::MatrixXd A0) {
  dd_Amatrix vRepMatrix = vRep->matrix;
  Eigen::MatrixXd M(vRep->rowsize, vRep->colsize - 1);
  for (unsigned int i = 0; i < vRep->rowsize; i++) {
    for (unsigned int j = 0; j < vRep->colsize - 1; j++) {
      M(i, j) = vRepMatrix[i][j + 1][0];
    }
  }
  dd_FreeMatrix(vRep);

  std::cout << "M:\n";
  std::cout << M << "\n\n";

  // move to origin
  for (Eigen::Index i = M.rows(); i-- > 0;) {
    M.row(i) -= M.row(0);
  }
  // rotate s.t. normal vector (=active subspace) is parallel to axis
  Eigen::MatrixXd Q;
  Eigen::HouseholderQR<Eigen::MatrixXd> qr(A0);
  Q = qr.householderQ();

  //  std::cout << "Q:\n" << Q << "\n\n";
  //  std::cout << "Q*A\n" << Q * A0 << "\n\n";

  M = (Q * M.transpose()).transpose();
  std::cout << "(Q * M^T)^T\n";
  std::cout << M << "\n---\n";
  return M;
}

double hPolytopeVolume(dd_MatrixPtr hRep) {
  typedef double NT;
  typedef Cartesian<double> Kernel;
  typedef boost::mt19937 RNGType;
  typedef HPolytope<Kernel::Point> Hpolytope;
  Hpolytope HP;
  // Create the polytope in form Ax<=b
  Hpolytope::MT LPA;
  Hpolytope::VT b;
  unsigned int dim = hRep->colsize - 1;
  dd_Amatrix hRepMatrix = hRep->matrix;

  LPA.resize(hRep->rowsize, hRep->colsize - 1);
  b.resize(hRep->rowsize);
  for (dd_rowrange i = 0; i < hRep->rowsize; i++) {
    b(i) = hRepMatrix[i][0][0];
    for (dd_colrange j = 1; j < hRep->colsize; j++) {
      LPA(i, j - 1) = hRepMatrix[i][j][0];
    }
  }
  dd_FreeMatrix(hRep);

  HP.init(dim, LPA, b);

  // Compute chebychev ball
  std::pair<Kernel::Point, NT> CheBall;
  CheBall = HP.ComputeInnerBall();

  // Parameter setup
  int n = HP.dimension();
  int walk_len = 10 + n / 10;
  int e = 1;  // goal error
  int rnum = std::pow(e, -2) * 400 * n * std::log(n);
  NT C = 2;
  NT ratio = 1.0 - 1.0 / (NT(n));
  int N = 500 * (static_cast<int>(C)) + (static_cast<int>(n * n / 2));
  int W = 4 * n * n + 500;

  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  RNGType rng(seed);
  boost::normal_distribution<> rdist(0, 1);
  boost::random::uniform_real_distribution<>(urdist);
  boost::random::uniform_real_distribution<> urdist1(-1, 1);

  vars<NT, RNGType> var1(rnum, n, walk_len, 1, 0, 1, 0, 0, 0, 0, rng, urdist, urdist1, -1.0, false,
                         false, false, false, false, false, true);
  // this consumes incredible amounts of memory (test for other cases)
  //  NT vol1 = volume(HP, var1, var1, CheBall);

  vars_g<NT, RNGType> var2(n, walk_len, N, W, 1, 0.2, CheBall.second, rng, C, 0.1, ratio, -1, false,
                           false, false, false, false, false, false, true);
  NT vol2 = volume_gaussian_annealing(HP, var2, var1, CheBall);

  return vol2;
}

int main() {
  dd_set_global_constants(); /* First, this must be called. */

  size_t originalDim = 2;
  Eigen::MatrixXd A0(2, 1);
  A0(0, 0) = 1.0 / sqrt(2.0);
  A0(1, 0) = 1.0 / sqrt(2.0);
  Eigen::VectorXd point(originalDim);
  point(0) = 1.0 / 3.0;
  point(1) = 0.0;
  dd_MatrixPtr hMatrix = createHPolytopeDescription(originalDim, A0, point);

  //  size_t originalDim = 3;
  //  Eigen::MatrixXd A0(originalDim, 1);
  //  A0(0, 0) = -1.0 / sqrt(3.0);
  //  A0(1, 0) = 1.0 / sqrt(3.0);
  //  A0(2, 0) = 1.0 / sqrt(3.0);
  //  Eigen::VectorXd point(originalDim);
  //  point(0) = 1;
  //  point(1) = 0;
  //  point(2) = 1;
  //  dd_MatrixPtr hMatrix = createHPolytopeDescription(originalDim, A0, point);

  dd_ErrorType err;
  dd_PolyhedraPtr poly = dd_DDMatrix2Poly(hMatrix, &err);
  dd_FreeMatrix(hMatrix);
  dd_Amatrix A = poly->A;
  std::cout << "H-Representation:\n";
  for (dd_rowrange m = 0; m < poly->m; m++) {
    for (dd_colrange d = 0; d < poly->d; d++) {
      std::cout << A[m][d][0] << " ";
    }
    std::cout << "\n";
  }
  dd_MatrixPtr vRep = dd_CopyGenerators(poly);
  dd_FreePolyhedra(poly);
  //  std::cout << "V-Representation:\n";
  //  dd_Amatrix vRepMatrix = vRep->matrix;
  //  for (size_t i = 0; i < vRep->rowsize; i++) {
  //    for (size_t j = 0; j < vRep->colsize; j++) {
  //      std::cout << vRepMatrix[i][j][0] << " ";
  //    }
  //    std::cout << "\n";
  //  }

  std::cout << "--------------\n";

  Eigen::MatrixXd M = affineHull(vRep, A0);

  // V -> H
  // delete axis parallel coordinate which is now obsolete (axes are enumbered 0,1,..)
  size_t obsoleteCoordinate = 0;
  std::cout << "removing " << obsoleteCoordinate << "'th coordinate.\n";
  dd_MatrixPtr vMatrix = dd_CreateMatrix(M.rows(), M.cols());
  for (size_t i = 0; i < M.rows(); i++) {
    vMatrix->matrix[i][0][0] = 1.0;
    for (size_t j = 0; j < M.cols(); j++) {
      if (j < obsoleteCoordinate) {
        vMatrix->matrix[i][j + 1][0] = M(i, j);
      } else if (j > obsoleteCoordinate) {
        vMatrix->matrix[i][j][0] = M(i, j);
      }
    }
  }

  std::cout << "vMatrix\n";
  for (size_t i = 0; i < M.rows(); i++) {
    for (size_t j = 0; j < M.cols(); j++) {
      std::cout << vMatrix->matrix[i][j][0] << " ";
    }
    std::cout << "\n";
  }

  vMatrix->numbtype = dd_Real;
  vMatrix->representation = dd_Generator;
  dd_PolyhedraPtr vPoly = dd_DDMatrix2Poly(vMatrix, &err);
  dd_FreeMatrix(vMatrix);

  // Eigentlich unnötig. Volesti kann auch direkt aus V-Darstellung das Volumen berechnen!
  std::cout << "H-Representation:\n";
  dd_MatrixPtr hRep = dd_CopyInequalities(vPoly);
  dd_FreePolyhedra(vPoly);
  dd_Amatrix hRepMatrix = hRep->matrix;
  for (dd_rowrange i = 0; i < hRep->rowsize; i++) {
    for (dd_colrange j = 0; j < hRep->colsize; j++) {
      std::cout << hRepMatrix[i][j][0] << " ";
    }
    std::cout << "\n";
  }

  std::cout << "--------------\n";
  double vol2 = hPolytopeVolume(hRep);
  std::cout << "Computed volume (alg.2) = " << vol2 << std::endl;
  return 0;
}
