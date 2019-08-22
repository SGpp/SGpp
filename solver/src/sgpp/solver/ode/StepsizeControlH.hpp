// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef STEPSIZECONTROLH_HPP
#define STEPSIZECONTROLH_HPP

#include <sgpp/base/application/ScreenOutput.hpp>
#include <sgpp/solver/ODESolver.hpp>
#include <sgpp/solver/ode/StepsizeControl.hpp>

#include <sgpp/globaldef.hpp>

#include <string>

namespace sgpp {
namespace solver {

/**
 * This class implements a step size control using Crank-Nicolson with different step sizes
 * for solving ordinary partial equations
 *
 */
class StepsizeControlH : public StepsizeControl {
 private:
  void predictor(SLESolver& LinearSystemSolver,
                 sgpp::solver::OperationParabolicPDESolverSystem& System, double tmp_timestepsize,
                 sgpp::base::DataVector& dv, sgpp::base::DataVector& corr,
                 sgpp::base::DataVector* rhs);
  void corrector(SLESolver& LinearSystemSolver,
                 sgpp::solver::OperationParabolicPDESolverSystem& System, double tmp_timestepsize,
                 sgpp::base::DataVector& dv, sgpp::base::DataVector* rhs);

  // double twoNorm(sgpp::base::DataVector &dv1, sgpp::base::DataVector &dv2);

  double nextTimestep(double tmp_timestepsize, double tmp_timestepsize_old, double norm,
                       double epsilon);

  std::string _odesolver;

 public:
  /**
   * Std-Constructer
   *
   * @param odesolver the mode of the euler that should be executed, must be ExEul or ImEul
   * @param imax number of maximum executed iterations
   * @param timestepSize the size of one timestep
   * @param eps the epsilon for the step size control
   * @param screen possible pointer to a sgpp::base::ScreenOutput object
   * @param gamma damping factor
   */
  StepsizeControlH(std::string odesolver, size_t imax, double timestepSize, double eps,
                   sgpp::base::ScreenOutput* screen = nullptr, double gamma = 0.9);

  /**
   * Std-Destructor
   */
  virtual ~StepsizeControlH();
};

}  // namespace solver
}  // namespace sgpp

#endif /* STEPSIZECONTROLH_HPP */
