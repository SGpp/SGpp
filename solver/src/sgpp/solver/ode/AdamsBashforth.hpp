// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef ADAMSBASHFORTH_HPP
#define ADAMSBASHFORTH_HPP

#include <sgpp/base/application/ScreenOutput.hpp>
#include <sgpp/solver/ODESolver.hpp>
#include <sgpp/globaldef.hpp>

#include <string>

namespace sgpp {
namespace solver {

/**
 * This class implements the explicit Adams-Bashforth method
 * for solving ordinary partial equations
 *
 */
class AdamsBashforth : public ODESolver {
 private:
  /// Pointer to sgpp::base::ScreenOutput object
  sgpp::base::ScreenOutput* myScreen;

 public:
  /**
   * Std-Constructer
   *
   * @param imax number of maximum executed iterations
   * @param timestepSize the size of one timestep
   * @param screen possible pointer to a sgpp::base::ScreenOutput object
   */
  AdamsBashforth(size_t imax, double timestepSize, sgpp::base::ScreenOutput* screen = nullptr);

  /**
   * Std-Destructor
   */
  virtual ~AdamsBashforth();

  virtual void solve(SLESolver& LinearSystemSolver,
                     sgpp::solver::OperationParabolicPDESolverSystem& System,
                     bool bIdentifyLastStep = false, bool verbose = false);
};

}  // namespace solver
}  // namespace sgpp

#endif /* ADAMSBASHFORTH_HPP */
