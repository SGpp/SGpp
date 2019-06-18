// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/datadriven/algorithm/DBMatDecompMatrixSolver.hpp>
#include <sgpp/datadriven/configuration/ParallelConfiguration.hpp>
#include <sgpp/datadriven/scalapack/DataMatrixDistributed.hpp>
#include <sgpp/datadriven/scalapack/DataVectorDistributed.hpp>

namespace sgpp {
namespace datadriven {

/**
 * This class solves an (lhs + lambda*I) * alpha = b system of linear equations
 * after the offline and online phases are done
 */

class DBMatDMS_SMW : public DBMatDecompMatrixSolver {
 public:
  /**
   * Constructor
   */
  DBMatDMS_SMW() = default;

  /**
   * Solves the system after the corresponding offline and online objects are
   * done with decomposing and adaptivity, resp.
   * The computation done: alpha = A_inv*b + B*b
   *
   * @param A_inv Inverse of the original lhs matrix
   * @param B     Storage of the online objects refined/coarsened points
   * @param b     The right side of the system
   * @param alpha The solution vector of the system, computed values go there
   */
  void solve(sgpp::base::DataMatrix& A_inv, sgpp::base::DataMatrix& B, sgpp::base::DataVector& b,
             sgpp::base::DataVector& alpha);

  /**
   * Parallel/Distributed version of solve.
   *
   * @param A_inv Inverse of a tridiagonal matrix
   * @param B     Storage of the online objects refined/coarsened points
   * @param b     The right side of the system
   * @param alpha The solution vector of the system, computed values go there
   */
  void solveParallel(DataMatrixDistributed& A_inv, DataMatrixDistributed& B,
                     DataVectorDistributed& b, DataVectorDistributed& alpha);
};

}  // namespace datadriven
}  // namespace sgpp
