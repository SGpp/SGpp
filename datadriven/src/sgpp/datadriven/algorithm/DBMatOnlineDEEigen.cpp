// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifdef USE_GSL

#include <sgpp/datadriven/algorithm/DBMatOnlineDEEigen.hpp>

#include <sgpp/datadriven/algorithm/DBMatDMSEigen.hpp>

namespace sgpp {
namespace datadriven {
DBMatOnlineDEEigen::DBMatOnlineDEEigen(DBMatOffline& offline, Grid& grid, double lambda,
    double beta)
    : DBMatOnlineDE{offline, grid, lambda, beta} {}

void DBMatOnlineDEEigen::solveSLE(DataVector& alpha, DataVector& b, Grid& grid,
    DensityEstimationConfiguration& densityEstimationConfig, bool do_cv) {
  DataMatrix& lhsMatrix = offlineObject.getDecomposedMatrix();

  // Solve the system:
  alpha.resizeZero(lhsMatrix.getNcols());

  size_t n = lhsMatrix.getNcols();
  DataVector e(n);
  lhsMatrix.getRow(n, e);
  DBMatDMSEigen esolver;

  esolver.solve(lhsMatrix, e, alpha, b, lambda);
}

} /* namespace datadriven */
} /* namespace sgpp */

#endif /*USE_GSL*/
