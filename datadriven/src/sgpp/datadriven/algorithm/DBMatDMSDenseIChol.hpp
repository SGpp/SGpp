// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

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
 * incomplete cholesky factorization on a dense matrix
 */
class DBMatDMSDenseIChol : public DBMatDMSChol {
 public:
  DBMatDMSDenseIChol(
      const sgpp::datadriven::DensityEstimationConfiguration& densityEstimationConfig,
      Grid& grid, double lambda, bool doCV);

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
   * Perform backward substitution solving the triangular system $A alpha = y$ with a parallel
   * Jaccobi solver.
   * @param decompMatrix Triangular matrix
   * @param y right hand side obtained by forward substitution
   * @param alpha the vector of unknowns we solve for
   */
  void choleskyBackwardSolve(const DataMatrix& decompMatrix, const DataVector& y,
                             DataVector& alpha) const override;

  /**
   * Perform forward substitution solving the triangular system $L y = b$ with a parallel
   * Jaccobi solver.
   * @param decompMatrix Triangular matrix
   * @param b right hand side of our initial system matrix we solve for
   * @param y the vector of unknowns we solve for
   */
  void choleskyForwardSolve(const DataMatrix& decompMatrix, const DataVector& b,
                            DataVector& y) const override;

 private:
  /**
   * update the mutable proxy object to avoid costly copy operations when modifying lambda.
   */
  void updateProxyMatrixLambda(double lambda_up) const;

  /**
   * Configuration of the matrix decomposition. In there we find the parameters
   * to configure the amount of sweeps for the parallel algorithms
   */
  const sgpp::datadriven::DensityEstimationConfiguration& densityEstimationConfig;

  /**
   * proxy object to avoid costly copy operations when modifying lambda.
   */
  mutable DataMatrix proxyMatrix;
};

} /* namespace datadriven */
} /* namespace sgpp */
