/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * DBMatDMSDenseIChol.cpp
 *
 *  Created on: Apr 16, 2017
 *      Author: Michael Lettrich
 */

#include <sgpp/datadriven/algorithm/DBMatDMSDenseIChol.hpp>

#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/operation/hash/OperationMatrix.hpp>
#include <sgpp/pde/operation/PdeOpFactory.hpp>

namespace sgpp {
namespace datadriven {
using sgpp::base::DataMatrix;
using sgpp::base::Grid;
using sgpp::base::OperationMatrix;

DBMatDMSDenseIChol::DBMatDMSDenseIChol(
    const sgpp::datadriven::DensityEstimationConfiguration& densityEstimationConfig, Grid& grid,
    double lambda, bool doCV)
    : DBMatDMSChol{}, densityEstimationConfig{densityEstimationConfig}, proxyMatrix{} {
  // initialize proxy matrix if we do cv
  if (doCV) {
    auto size = grid.getStorage().getSize();
    proxyMatrix.resizeQuadratic(size);
    std::unique_ptr<OperationMatrix> op(
        sgpp::op_factory::createOperationLTwoDotExplicit(&proxyMatrix, grid));

    // set regularization parameter
    updateProxyMatrixLambda(lambda);
  }
}

void DBMatDMSDenseIChol::choleskyUpdateLambda(sgpp::base::DataMatrix& decompMatrix,
                                              double lambdaUpdate) const {
  updateProxyMatrixLambda(lambdaUpdate);
  DBMatOfflineDenseIChol::ichol(proxyMatrix, decompMatrix,
                                densityEstimationConfig.iCholSweepsUpdateLambda_);
}

void DBMatDMSDenseIChol::choleskyBackwardSolve(const sgpp::base::DataMatrix& decompMatrix,
                                               const sgpp::base::DataVector& y,
                                               sgpp::base::DataVector& alpha) const {
  // cache efficient version of jaccobi based backward substitution
  DataVector tmpVec{alpha.getSize()};
  alpha.setAll(0.0);
  auto size = alpha.getSize();

#pragma omp parallel
  {
    for (auto sweep = 0u; sweep < densityEstimationConfig.iCholSweepsSolver_; sweep++) {
      tmpVec.setAll(0.0);
#pragma omp for schedule(guided) nowait
      for (auto i = 0u; i < size; i++) {
#pragma omp simd
        for (auto j = 0u; j < i; j++) {
          tmpVec[j] += decompMatrix.get(i, j) * alpha.get(i);
        }
      }
#pragma omp parallel
      for (auto i = 0u; i < size; i++) {
        alpha.set(i, 1.0 / decompMatrix.get(i, i) * (y.get(i) - tmpVec.get(i)));
      }
    }
  }
}

void DBMatDMSDenseIChol::choleskyForwardSolve(const sgpp::base::DataMatrix& decompMatrix,
                                              const sgpp::base::DataVector& b,
                                              sgpp::base::DataVector& y) const {
  // initial guess for y
  y.setAll(0.0);

  auto size = y.getSize();

#pragma omp parallel
  {
    for (auto sweep = 0u; sweep < densityEstimationConfig.iCholSweepsSolver_; sweep++) {
#pragma omp for schedule(guided) nowait
      for (auto i = 0u; i < size; i++) {
        auto tmp = 0.0;
        for (auto j = 0u; j < i; j++) {
          tmp += decompMatrix.get(i, j) * y.get(j);
        }
        y.set(i, 1.0 / decompMatrix.get(i, i) * (b.get(i) - tmp));
      }
    }
  }
}

void DBMatDMSDenseIChol::updateProxyMatrixLambda(double lambdaUpdate) const {
  auto size = proxyMatrix.getNrows();
  for (auto i = 0u; i < size; i++) {
    auto value = proxyMatrix.get(i, i) + lambdaUpdate;
    proxyMatrix.set(i, i, value);
  }
}

} /* namespace datadriven */
} /* namespace sgpp */
