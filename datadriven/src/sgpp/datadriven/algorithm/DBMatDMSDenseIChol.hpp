/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * DBMatDMSDenseIChol.hpp
 *
 *  Created on: Apr 16, 2017
 *      Author: Michael Lettrich
 */

#pragma once
#include <sgpp/datadriven/algorithm/DBMatDMSChol.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineDenseIChol.hpp>

#include <sgpp/base/grid/Grid.hpp>

namespace sgpp {
namespace datadriven {

using sgpp::base::Grid;
using sgpp::base::DataMatrix;
using sgpp::base::DataVector;

/**
 * Solve the system of equations with a LL'-decomposed matrix where LL' is created by an iterative,
 * incomplete
 * cholesky factorization on a dense matrix
 */
class DBMatDMSDenseIChol : public DBMatDMSChol {
 public:
  DBMatDMSDenseIChol(const DBMatOfflineIcholParameters& params, Grid& grid, double lambda,
                     bool doCV);

 protected:
  /**
   * Update the regularization factor of the decomposition. This is a very costly operation as the
   * iterative algorithm will need to recalculate the entire decomposition. However the current
   * solution is used as an initial guess, reducing the required amount of sweeps.
   * @param decompMatrix
   * @param lambdaUp
   */
  void choleskyUpdateLambda(DataMatrix& decompMatrix, double lambdaUp) const override;

  /**
   * Perform backward substit
   */
  void choleskyBackwardSolve(const DataMatrix& decompMatrix, const DataVector& y,
                             DataVector& alpha) const override;

  /**
   * TODO(lettrich) : write documentation
   */
  void choleskyForwardSolve(const DataMatrix& decompMatrix, const DataVector& b,
                            DataVector& y) const override;

 private:
  /**
   * TODO(lettrich) : write documentation
   */
  void updateProxyMatrixLambda(double lambda_up) const;

  /**
   * TODO(lettrich) : write documentation
   */
  DBMatOfflineIcholParameters params;

  /**
   * TODO(lettrich) : write documentation
   */
  mutable DataMatrix proxyMatrix;
};

} /* namespace datadriven */
} /* namespace sgpp */
