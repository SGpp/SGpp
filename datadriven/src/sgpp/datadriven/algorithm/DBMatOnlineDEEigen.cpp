/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * DBMatOnlineDEEigen.cpp
 *
 *  Created on: Apr 8, 2017
 *      Author: Michael Lettrich
 */

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

  if (canCV && do_cv) {
    double best_crit = 0;
    double cur_lambda;
    for (int i = 0; i < lambdaStep; i++) {
      cur_lambda = lambdaStart + i * (lambdaEnd - lambdaStart) / (lambdaStep - 1);
      if (cvLogscale) cur_lambda = exp(cur_lambda);
      esolver.solve(lhsMatrix, e, alpha, b, cur_lambda);
      // double crit = computeL2Error();
      double crit = resDensity(alpha, grid);
      // std::cout << "cur_lambda: " << cur_lambda << ", crit: " << crit <<
      // std::endl;
      if (i == 0 || crit < best_crit) {
        best_crit = crit;
        lambda = cur_lambda;
      }
    }
  }
  esolver.solve(lhsMatrix, e, alpha, b, lambda);
}

} /* namespace datadriven */
} /* namespace sgpp */

#endif /*USE_GSL*/
