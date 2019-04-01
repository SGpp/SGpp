// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/datadriven/algorithm/DBMatDecompMatrixSolver.hpp>
#include <sgpp/datadriven/scalapack/DataMatrixDistributed.hpp>
#include <sgpp/datadriven/scalapack/DataVectorDistributed.hpp>

namespace sgpp {
namespace datadriven {

/**
 * Class to solve the system of equations with a LL'-decomposed matrix
 */

class DBMatDMSChol : public DBMatDecompMatrixSolver {
 public:
  /**
   * Default constructor
   */
  DBMatDMSChol() = default;

  /**
   * Solves a system of equations
   *
   * @param decompMatrix the LL' lower triangular cholesky factor
   * @param alpha the vector of unknowns (the result is stored there)
   * @param b the right hand vector of the equation system
   * @param lambda_old the current regularization paramter
   * @param lambda_new the new regularization paramter (e.g. if cross-validation
   * is applied)
   */
  virtual void solve(sgpp::base::DataMatrix& decompMatrix, sgpp::base::DataVector& alpha,
                     const sgpp::base::DataVector& b, double lambda_old, double lambda_new) const;

  /**
   * Parallel (distributed) version of solve.
   * @param decompMatrix the LL' lower triangular cholesky factor
   * @param x input: the right hand vector of the equation system, output: the vector of unknowns
   * (the result is stored there)
   * @param lambda_old the current regularization paramter
   * @param lambda_new the new regularization paramter (e.g. if cross-validation
   * is applied)
   */
  void solveParallel(DataMatrixDistributed& decompMatrix, DataVectorDistributed& alpha,
                     double lambda_old, double lambda_new) const;

  /**
   * Performe a rank one cholesky update
   *
   * @param decompMatrix the LL' lower triangular cholesky factor
   * @param update the vector representing the rank one matrix (xx')
   * @param do_cv indicating if updates are used for cross valdiation
   * (using special structure of update vectors to save runtime)
   */
  void choleskyUpdate(sgpp::base::DataMatrix& decompMatrix, const sgpp::base::DataVector& update,
                      bool do_cv = false) const;

  /**
   * Performe a rank one cholesky downdate
   *
   * @param decompMatrix the LL' lower triangular cholesky factor
   * @param downdate the vector representing the rank one matrix (xx')
   * @param do_cv indicating if updates are used for cross valdiation
   * (using special structure of update vectors to save runtime)
   */
  void choleskyDowndate(sgpp::base::DataMatrix& decompMatrix,
                        const sgpp::base::DataVector& downdate, bool do_cv = false) const;

 protected:
  /**
   * Update the decomposition if the regularization parameter changes. This may be more expensive
   * then recalculating the decomposition.
   * @param decompMatrix decomposed matrix to be modified
   * @param lambdaUpdate the value by which the regularization parameter modifies the diagonal.
   */
  virtual void choleskyUpdateLambda(sgpp::base::DataMatrix& decompMatrix,
                                    double lambdaUpdate) const;

  /**
   * Perform Backward substitution solving the triangular system $A alpha = y$
   * @param decompMatrix Triangular matrix
   * @param y right hand side obtained by forward substitution
   * @param alpha the vector of unknowns we solve for
   */
  virtual void choleskyBackwardSolve(const sgpp::base::DataMatrix& decompMatrix,
                                     const sgpp::base::DataVector& y,
                                     sgpp::base::DataVector& alpha) const;

  /**
   * Perform forward substitution solving the triangular system $L y = b$
   * @param decompMatrix Triangular matrix
   * @param b right hand side of our initial system matrix we solve for
   * @param y the vector of unknowns we solve for
   */
  virtual void choleskyForwardSolve(const sgpp::base::DataMatrix& decompMatrix,
                                    const sgpp::base::DataVector& b,
                                    sgpp::base::DataVector& y) const;
};

}  // namespace datadriven
}  // namespace sgpp
