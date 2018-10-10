// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/optimization/activeSubspaces/ASMatrix.hpp>

namespace sgpp {
namespace optimization {

void ASMatrix::evDecompositionForSymmetricMatrices() {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(C);
  if (eigensolver.info() != Eigen::Success) abort();
  this->eigenvalues = eigensolver.eigenvalues();

  this->W = eigensolver.eigenvectors();
}

void ASMatrix::setMatrix(Eigen::MatrixXd newC) {
  if ((newC.cols() == static_cast<unsigned int>(numDim)) &&
      (newC.rows() == static_cast<unsigned int>(numDim))) {
    C = newC;
  } else {
    std::cout << "ASMatrix: matrix size does not match objective function" << std::endl;
  }
}

}  // namespace optimization
}  // namespace sgpp
