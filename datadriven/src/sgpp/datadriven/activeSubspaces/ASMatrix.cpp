// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/activeSubspaces/ASMatrix.hpp>

namespace sgpp {
namespace datadriven {

void ASMatrix::evDecompositionForSymmetricMatrices() {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(C);
  if (eigensolver.info() != Eigen::Success) {
    std::cout << "ASMatrix: eigenvalue decomposition failed\n";
    abort();
  }
  this->eigenvalues = eigensolver.eigenvalues();
  W = eigensolver.eigenvectors();

  // We don't want all entries of an eigenvecctor  Wi to be negative because this leads to problems
  // with leftBound1D and rightBound1D in ASResponseSurfaceNakBspline.
  // Heuristic: if one arbitrary entry is negative simply turn around Wi => not all entries are
  // negative
  // ToDo(rehmemk) Is this safe in higher dimensions / are there any problems with the bounds?
  for (unsigned int i = 0; i < W.cols(); i++) {
    if (W(0, i) < 0) {
      W.col(i) *= -1;
      break;
    }
  }
}

}  // namespace datadriven
}  // namespace sgpp
