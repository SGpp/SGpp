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

#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/algorithm/DBMatDMSDenseIChol.hpp>

#include <iomanip>
#include <list>
#include <vector>

namespace sgpp {
namespace datadriven {

DBMatOnlineDEChol::DBMatOnlineDEChol(DBMatOffline& offline, double beta)
    : DBMatOnlineDE{offline, beta} {}

void DBMatOnlineDEChol::solveSLE(DataVector& b, bool do_cv) {
  DataMatrix& lhsMatrix = offlineObject.getDecomposedMatrix();
  alpha = DataVector(lhsMatrix.getNcols());

  auto cholsolver = std::unique_ptr<DBMatDMSChol>{buildCholSolver(offlineObject, do_cv)};

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
      cholsolver->solve(lhsMatrix, alpha, b, old_lambda, cur_lambda);
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
  cholsolver->solve(lhsMatrix, alpha, b, old_lambda, lambda);

  //  DBMatDMSChol myCholSolver;
  //  DataVector myAlpha{alpha.getSize()};
  //  myCholSolver.solve(lhsMatrix, myAlpha, b, old_lambda, lambda);
  //
  //  myAlpha.sub(alpha);
  //  myAlpha.abs();
  //  myAlpha.sqr();
  //  auto res = sqrt(myAlpha.sum());
  //  std::cout << "solving with " << offlineObject.getDensityEstimationConfig().icholSweepsSolver_
  //            << " sweeps results in error: " << std::scientific << std::setprecision(10) << res
  //            << "\n";
}

DBMatDMSChol* DBMatOnlineDEChol::buildCholSolver(DBMatOffline& offlineObject, bool doCV) const {
  // const cast is OK here, since we access the config read only.
  switch (offlineObject.getDensityEstimationConfig().decomposition_) {
    case (MatrixDecompositionType::Chol):
      return new DBMatDMSChol();
      break;
    case (MatrixDecompositionType::DenseIchol):
      return new DBMatDMSDenseIChol(offlineObject.getDensityEstimationConfig(),
                                    offlineObject.getGrid(),
                                    offlineObject.getRegularizationConfig().lambda_,
                                    doCV);
      break;
    case (MatrixDecompositionType::LU):
    case (MatrixDecompositionType::Eigen):
    case (MatrixDecompositionType::OrthoAdapt):
    default:
      throw sgpp::base::algorithm_exception{"Only Cholesky based solvers can use this Solver"};
  }
}

std::vector<size_t> DBMatOnlineDEChol::updateSystemMatrixDecomposition(
    size_t numAddedGridPoints,
    std::list<size_t> deletedGridPointIndices,
    double lambda)  {
  DBMatOffline* offlineObject = &getOfflineObject();
  dynamic_cast<DBMatOfflineChol *>(offlineObject)
      ->choleskyModification(numAddedGridPoints, deletedGridPointIndices, lambda);
      std::vector<size_t> return_vector = {};
      return return_vector;
}


} /* namespace datadriven */
} /* namespace sgpp */
