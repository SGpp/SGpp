// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifdef USE_GSL

#ifndef DBMatDMSChol_HPP_
#define DBMatDMSChol_HPP_

#include <sgpp/datadriven/algorithm/DBMatDecompMatrixSolver.hpp>

/**
 * Class to solve the system of equations with a LU-decomposed matrix
 */

class DBMatDMSChol: public DBMatDecompMatrixSolver {
 public:
  /**
   * (Empty) constructor
   */
  DBMatDMSChol();

  /**
   * (Empty) destructor
   */
  virtual ~DBMatDMSChol();

  /**
   * Solves a system of equations
   *
   * @param DecompMatrix the LL' lower triangular cholesky factor
   * @param alpha the vector of unknowns (the result is stored there)
   * @param b the right hand vector of the equation system
   * @param lambda_old the current regularization paramter 
   * @param lambda_new the new regularization paramter (e.g. if cross-validation is applied)
   */
  void solve(sgpp::base::DataMatrix& DecompMatrix, sgpp::base::DataVector& alpha,
	     sgpp::base::DataVector& b, double lambda_old, double lambda_new);

  /**
   * Performe a rank one cholesky update 
   *
   * @param DecompMatrix the LL' lower triangular cholesky factor
   * @param update the vector representing the rank one matrix (xx')
   * @param do_cv indicating if updates are used for cross valdiation 
   * (using special structure of update vectors to save runtime)
   */
  void choleskyUpdate(sgpp::base::DataMatrix& DecompMatrix, 
                      sgpp::base::DataVector* update, bool do_cv = false);

  /**
   * Performe a rank one cholesky downdate 
   *
   * @param DecompMatrix the LL' lower triangular cholesky factor
   * @param downdate, the vector representing the rank one matrix (xx')
   * @param do_cv indicating if updates are used for cross valdiation 
   * (using special structure of update vectors to save runtime)
   */
  void choleskyDowndate(sgpp::base::DataMatrix& DecompMatrix,
                        sgpp::base::DataVector* downdate, bool do_cv = false);
};

#endif /* DBMatDMSChol_HPP_ */

#endif /* USE_GSL */


