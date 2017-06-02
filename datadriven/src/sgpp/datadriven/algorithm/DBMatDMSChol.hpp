// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifdef USE_GSL

#ifndef DBMatDMSChol_HPP_
#define DBMatDMSChol_HPP_

#include <sgpp/datadriven/algorithm/DBMatDecompMatrixSolver.hpp>

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
   * @param downdate, the vector representing the rank one matrix (xx')
   * @param do_cv indicating if updates are used for cross valdiation
   * (using special structure of update vectors to save runtime)
   */
  void choleskyDowndate(sgpp::base::DataMatrix& decompMatrix,
                        const sgpp::base::DataVector& downdate, bool do_cv = false) const;

 protected:
  /**
   * TODO(lettrich) : write documentation
   */
  virtual void choleskyUpdateLambda(sgpp::base::DataMatrix& decompMatrix, double lambda_up) const;

  /**
   * TODO(lettrich) : write documentation
   */
  virtual void choleskyBackwardSolve(const sgpp::base::DataMatrix& decompMatrix,
                                     const sgpp::base::DataVector& y,
                                     sgpp::base::DataVector& alpha) const;

  /**
   * TODO(lettrich) : write documentation
   */
  virtual void choleskyForwardSolve(const sgpp::base::DataMatrix& decompMatrix,
                                    const sgpp::base::DataVector& b,
                                    sgpp::base::DataVector& y) const;
};

}  // namespace datadriven
}  // namespace sgpp

#endif /* DBMatDMSChol_HPP_ */

#endif /* USE_GSL */
