// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_SLE_SOLVER_ARMADILLO_HPP
#define SGPP_OPTIMIZATION_SLE_SOLVER_ARMADILLO_HPP

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/sle/solver/SLESolver.hpp>

#include <vector>

namespace SGPP {
namespace optimization {
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
  bool solve(SLE& system, base::DataVector& b, base::DataVector& x) const override;

  /**
   * @param       system  system to be solved
   * @param       B       matrix of right-hand sides
   * @param[out]  X       matrix of solutions to the systems
   * @return              whether all went well
   *                      (false if errors occurred)
   */
  bool solve(SLE& system, base::DataMatrix& B, base::DataMatrix& X) const override;
};
}  // namespace sle_solver
}  // namespace optimization
}  // namespace SGPP

#endif /* SGPP_OPTIMIZATION_SLE_SOLVER_ARMADILLO_HPP */
