// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/tools/sle/solver/SLESolver.hpp>
#include <sgpp/globaldef.hpp>

#include <vector>

namespace sgpp {
namespace base {
namespace sle_solver {

/**
 * Automatic choice of external linear solver.
 */
class Auto : public SLESolver {
 public:
  /// maximal matrix dimension to allow use of full solvers
  static const size_t MAX_DIM_FOR_FULL = 30000;
  /// maximal matrix dimension to prefer GaussianElimination over BiCGStab
  static const size_t MAX_DIM_FOR_GAUSSIAN = 200;
  /// maximal ratio of non-zero entries for sparse solvers
  static constexpr double MAX_NNZ_RATIO_FOR_SPARSE = 0.1;
  /// ratio of the rows (e.g. 0.1 = 10%) to use for sparsity estimation
  static constexpr double ESTIMATE_NNZ_ROWS_SAMPLE_SIZE = 0.05;

  /**
   * Destructor.
   */
  ~Auto() override;

  /**
   * @param       system  system to be solved
   * @param       b       right-hand side
   * @param[out]  x       solution to the system
   * @return              whether all went well
   *                      (false if errors occurred)
   */
  bool solve(SLE& system, DataVector& b, DataVector& x) const override;

  /**
   * @param       system  system to be solved
   * @param       B       matrix of right-hand sides
   * @param[out]  X       matrix of solutions to the systems
   * @return              whether all went well
   *                      (false if errors occurred)
   */
  bool solve(SLE& system, DataMatrix& B, DataMatrix& X) const override;
};
}  // namespace sle_solver
}  // namespace base
}  // namespace sgpp
