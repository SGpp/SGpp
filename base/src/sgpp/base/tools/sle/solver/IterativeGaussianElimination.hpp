// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/tools/sle/solver/SLESolver.hpp>

#include <cstddef>
#include <stdexcept>
#include <vector>

namespace sgpp {
namespace base {
namespace sle_solver {

/**
 * @brief An incremental linear system solver based on Gaussian elimination.
 *
 * This solver computes a fully-pivoted LU decomposition of the system matrix.
 * For subsequent solves of systems that have been extended with new rows and columns,
 * it efficiently updates the existing decomposition rather than recomputing it from scratch.
 *
 * @warning The solver has no way of measuring whether numerical errors accumulate. However, if a
 * singularity is detected (pivot < tolerance) during an update, the solver automatically
 * performs a complete recalculation from scratch to ensure stability.
 */
class IterativeGaussianElimination : public SLESolver {
 public:
  /// Tolerance to detect near-zero pivots, which could indicate matrix degeneration.
  constexpr static const double DEGENERATION_TOLERANCE = 1e-9;

  /**
   * @brief Default constructor.
   */
  IterativeGaussianElimination();

  /**
   * @brief Destructor.
   */
  ~IterativeGaussianElimination() override;

  /**
   * This class maintains state for incremental updates, so the const version of solve() is not
   * supported. Calling this method will throw a std::runtime_error. Please use iterativeSolve()
   * instead.
   */
  bool solve(SLE& system, DataVector& b, DataVector& x) const override {
    throw std::runtime_error(
        "IterativeGaussianElimination::solve is not supported. Use iterativeSolve().");
  }

  /**
   * Solves the linear system Ax = b using an incremental LU decomposition.
   *
   * If the system matrix A has grown since the last call, the solver will update
   * its internal LU decomposition. If the update fails due to numerical instability
   * (singularity detected), a full decomposition is computed from scratch.
   *
   * @warning This implementation assumes that the existing (old) part of the
   * system matrix does not change between calls. If it does change, or if the
   * system shrinks, a full re-solve is automatically triggered.
   *
   * @param      system  The linear system to be solved (provides the matrix A).
   * @param      b       The right-hand side vector.
   * @param[out] x       The solution vector.
   * @return             True if the solution was successful, false otherwise (e.g., matrix is
   * singular).
   */
  bool iterativeSolve(SLE& system, DataVector& b, DataVector& x);

 protected:
  /// A copy of the system matrix from the last solve.
  DataMatrix A;
  /// The in-place LU decomposition of the matrix A.
  DataMatrix LU;
  /// Pivot indices for rows.
  std::vector<size_t> pivotRow;
  /// Pivot indices for columns (used for full pivoting).
  std::vector<size_t> pivotCol;
};
}  // namespace sle_solver
}  // namespace base
}  // namespace sgpp
