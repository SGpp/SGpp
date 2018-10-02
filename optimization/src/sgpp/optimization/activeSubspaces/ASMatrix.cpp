// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/optimization/activeSubspaces/ASMatrix.hpp>

namespace sgpp {
namespace optimization {

void ASMatrix::evDecomposition() {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(this->C);
  if (eigensolver.info() != Eigen::Success) abort();
  this->eigenvalues = eigensolver.eigenvalues();
  this->W = eigensolver.eigenvectors();
}

}  // namespace optimization
}  // namespace sgpp
