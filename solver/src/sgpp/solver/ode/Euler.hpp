// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef EULER_HPP
#define EULER_HPP

#include <sgpp/base/application/ScreenOutput.hpp>
#include <sgpp/solver/ODESolver.hpp>
#include <sgpp/globaldef.hpp>

#include <string>

namespace sgpp {
namespace solver {

/**
 * This class implements the explicit and implicit Euler method
 * for solving ordinary partial equations
 *
 */
class Euler : public ODESolver {
 private:
  /// specifies if a grid evaluation should be execute in every time step
  bool bAnimation;
  /// specifies the type of euler that should be executed
  std::string ExMode;
  /// Pointer to sgpp::base::ScreenOutput object
  sgpp::base::ScreenOutput* myScreen;

 public:
  /**
   * Std-Constructer
   *
   * @param Mode the mode of the euler that should be executed, must be ExEul or ImEul
   * @param imax number of maximum executed iterations
   * @param timestepSize the size of one timestep
   * @param generateAnimation set this, if you want to create a grid evaluation in every time step,
   * in order to create an animation
   * @param screen possible pointer to a sgpp::base::ScreenOutput object
   */
  Euler(std::string Mode, size_t imax, double timestepSize, bool generateAnimation = false,
        sgpp::base::ScreenOutput* screen = nullptr);

  /**
   * Std-Destructor
   */
  virtual ~Euler();

  virtual void solve(SLESolver& LinearSystemSolver,
                     sgpp::solver::OperationParabolicPDESolverSystem& System,
                     bool bIdentifyLastStep = false, bool verbose = false);
};

}  // namespace solver
}  // namespace sgpp

#endif /* EULER_HPP */
