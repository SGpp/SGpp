// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/tools/sle/solver/SLESolver.hpp>

#include <cstddef>

namespace sgpp {
namespace base {
namespace sle_solver {

/**
 * Linear system solver implementing the direct Gaussian elimination.
 */
class GaussianElimination : public SLESolver {
 public:
  /**
   * Destructor.
   */
  ~GaussianElimination() override;

  /**
   * @param       system  system to be solved
   * @param       b       right-hand side
   * @param[out]  x       solution to the system
   * @return              whether all went well
   *                      (false if errors occurred)
   */
  bool solve(SLE& system, DataVector& b, DataVector& x) const override;
};
}  // namespace sle_solver
}  // namespace base
}  // namespace sgpp
