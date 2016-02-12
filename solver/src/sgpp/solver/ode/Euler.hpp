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

namespace SGPP {
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
  /// specifies the evaluation per dimension when a animation is created
  size_t evalsAnimation;
  /// specifies the type of euler that should be executed
  std::string ExMode;
  /// Pointer to SGPP::base::ScreenOutput object
  SGPP::base::ScreenOutput* myScreen;

 public:
  /**
   * Std-Constructer
   *
   * @param Mode the mode of the euler that should be executed, must be ExEul or ImEul
   * @param imax number of maximum executed iterations
   * @param timestepSize the size of one timestep
   * @param generateAnimation set this, if you want to create a grid evaluation in every time step,
   * in order to create an animation
   * @param numEvalsAnimation specifies the evaluation per dimension when a animation is created
   * @param screen possible pointer to a SGPP::base::ScreenOutput object
   */
  Euler(std::string Mode, size_t imax, float_t timestepSize, bool generateAnimation = false,
        size_t numEvalsAnimation = 20, SGPP::base::ScreenOutput* screen = NULL);

  /**
   * Std-Destructor
   */
  virtual ~Euler();

  virtual void solve(SLESolver& LinearSystemSolver,
                     SGPP::solver::OperationParabolicPDESolverSystem& System,
                     bool bIdentifyLastStep = false, bool verbose = false);
};

}  // namespace solver
}  // namespace SGPP

#endif /* EULER_HPP */
