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
 * Linear system solver using Armadillo (direct full solver).
 */
class Armadillo : public SLESolver {
 public:
  /**
   * Destructor.
   */
  ~Armadillo() override;

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
