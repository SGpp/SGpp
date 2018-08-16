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

DBMatOnlineDEChol::DBMatOnlineDEChol(DBMatOffline& offline, Grid& grid, double lambda, double beta)
    : DBMatOnlineDE{offline, grid, lambda, beta} {}

void DBMatOnlineDEChol::solveSLE(DataVector& alpha, DataVector& b, Grid& grid,
    DensityEstimationConfiguration& densityEstimationConfig, bool do_cv) {
  DataMatrix& lhsMatrix = offlineObject.getDecomposedMatrix();
  alpha.resizeZero(lhsMatrix.getNcols());

  auto cholsolver = std::unique_ptr<DBMatDMSChol>{buildCholSolver(offlineObject, grid,
      densityEstimationConfig, do_cv)};

  // Solve for density declaring coefficients alpha
  // std::cout << "lambda: " << lambda << std::endl;
  cholsolver->solve(lhsMatrix, alpha, b, lambda, lambda);

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

DBMatDMSChol* DBMatOnlineDEChol::buildCholSolver(DBMatOffline& offlineObject, Grid& grid,
    DensityEstimationConfiguration& densityEstimationConfig, bool doCV) const {
  // const cast is OK here, since we access the config read only.
  switch (offlineObject.getDecompositionType()) {
    case (MatrixDecompositionType::Chol):
      return new DBMatDMSChol();
      break;
    case (MatrixDecompositionType::DenseIchol):
      return new DBMatDMSDenseIChol(densityEstimationConfig,
                                    grid,
                                    lambda,
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
    DensityEstimationConfiguration& densityEstimationConfig,
    Grid& grid,
    size_t numAddedGridPoints,
    std::list<size_t> deletedGridPointIndices,
    double lambda)  {
  DBMatOffline* offlineObject = &getOfflineObject();
  dynamic_cast<DBMatOfflineChol *>(offlineObject)
      ->choleskyModification(grid, densityEstimationConfig, numAddedGridPoints,
          deletedGridPointIndices, lambda);
      std::vector<size_t> return_vector = {};
      return return_vector;
}


} /* namespace datadriven */
} /* namespace sgpp */
