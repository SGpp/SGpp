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

class DBMatDMSOrthoAdapt : public DBMatDecompMatrixSolver {
 public:
  /**
   * Constructor
   */
  DBMatDMSOrthoAdapt() = default;

  /**
   * Solves the system after the corresponding offline and online objects are
   * done with decomposing and adaptivity, resp.
   * The computation done: alpha = Q*T_inv*Q^t*b + B*b
   *
   * @param T_inv Inverse of a tridiagonal matrix
   * @param Q     Orthogonal matrix, part of hessenberg_decomp of the lhs matrix
   * @param B     Storage of the online objects refined/coarsened points
   * @param b     The right side of the system
   * @param alpha The solution vector of the system, computed values go there
   */
  void solve(sgpp::base::DataMatrix& T_inv, sgpp::base::DataMatrix& Q, sgpp::base::DataMatrix& B,
             sgpp::base::DataVector& b, sgpp::base::DataVector& alpha);

  /**
   * Parallel (distributed) version of solve.
   *
   * @param T_inv Inverse of a tridiagonal matrix
   * @param Q     Orthogonal matrix, part of hessenberg_decomp of the lhs matrix
   * @param B     Storage of the online objects refined/coarsened points
   * @param b     The right side of the system
   * @param alpha The solution vector of the system, computed values go there
   */
  void solveParallel(DataMatrixDistributed& T_inv, DataMatrixDistributed& Q,
                     DataMatrixDistributed& B, DataVectorDistributed& b,
                     DataVectorDistributed& alpha, std::shared_ptr<BlacsProcessGrid> processGrid,
                     const ParallelConfiguration& parallelConfig);
};

}  // namespace datadriven
}  // namespace sgpp
