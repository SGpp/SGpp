/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * DBMatOnlineDEChol.cpp
 *
 *  Created on: Apr 8, 2017
 *      Author: Michael Lettrich
 */

#include <sgpp/datadriven/algorithm/DBMatOnlineDEChol.hpp>

#include <sgpp/datadriven/algorithm/DBMatDMSChol.hpp>

namespace sgpp {
namespace datadriven {} /* namespace datadriven */
} /* namespace sgpp */

sgpp::datadriven::DBMatOnlineDEChol::DBMatOnlineDEChol(DBMatOffline& offline, double beta)
    : DBMatOnlineDE{offline, beta} {}

void sgpp::datadriven::DBMatOnlineDEChol::solveSLE(DataVector& b, bool do_cv) {
  DataMatrix& lhsMatrix = offlineObject.getDecomposedMatrix();
  alpha = DataVector(lhsMatrix.getNcols());
  DBMatDMSChol cholsolver;

  double old_lambda = lambda;
  // Perform cross-validation based on rank one up- and downdates
  // -> SHOULD NOT BE USED FOR LARGER GRID SETTINGS
  // ToDo: May be speed up by parallelization
  if (canCV && do_cv) {
    double best_crit = 0;
    double cur_lambda;
    for (int i = 0; i < lambdaStep; i++) {
      cur_lambda = lambdaStart + i * (lambdaEnd - lambdaStart) / (lambdaStep - 1);
      if (cvLogscale) cur_lambda = exp(cur_lambda);
      // std::cout << "Cur_lambda: " << cur_lambda << "  Old_lambda: " <<
      // old_lambda << std::endl;
      // Solve for density declaring coefficients alpha based on changed
      // lambda
      cholsolver.solve(lhsMatrix, alpha, b, old_lambda, cur_lambda);
      old_lambda = cur_lambda;
      double crit = resDensity(alpha);
      // std::cout << ", crit: " << crit << std::endl;
      if (i == 0 || crit < best_crit) {
        best_crit = crit;
        lambda = cur_lambda;
      }
    }
  }
  // Solve for density declaring coefficients alpha
  // std::cout << "lambda: " << lambda << std::endl;
  cholsolver.solve(lhsMatrix, alpha, b, old_lambda, lambda);
}
