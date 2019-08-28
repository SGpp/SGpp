// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef STEPSIZECONTROLEJ_HPP
#define STEPSIZECONTROLEJ_HPP

#include <sgpp/base/application/ScreenOutput.hpp>
#include <sgpp/solver/ODESolver.hpp>
#include <sgpp/solver/ode/StepsizeControl.hpp>

#include <sgpp/globaldef.hpp>

#include <string>

namespace sgpp {
namespace solver {

/**
 * This class implements a time step size control using 1D-Diffusion
 * for solving ordinary partial equations
 *
 * For solving the system of linear equations the
 * already implemented CG-method is used
 *
 */
class StepsizeControlEJ : public StepsizeControl {
 private:
  virtual void predictor(SLESolver& LinearSystemSolver,
                         sgpp::solver::OperationParabolicPDESolverSystem& System,
                         double tmp_timestepsize, sgpp::base::DataVector& dv,
                         sgpp::base::DataVector& corr, sgpp::base::DataVector* rhs);

  virtual void corrector(SLESolver& LinearSystemSolver,
                         sgpp::solver::OperationParabolicPDESolverSystem& System,
                         double tmp_timestepsize, sgpp::base::DataVector& dv,
                         sgpp::base::DataVector* rhs);

  virtual double nextTimestep(double tmp_timestepsize, double tmp_timestepsize_old, double norm,
                               double epsilon);

  virtual double norm(sgpp::solver::OperationParabolicPDESolverSystem& System,
                       sgpp::base::DataVector& dv1, sgpp::base::DataVector& dv2);
  std::string _odesolver;

 public:
  /**
   * Std-Constructer
   *
   * @param odesolver the selected ode solver
   * @param nTimesteps number of maximum executed iterations
   * @param timestepSize the size of one timestep
   * @param eps the epsilon for the stepsize control
   * @param sc
   * @param screen possible pointer to a sgpp::base::ScreenOutput object
   * @param gamma used damping factor, default is 0.5
   */
  StepsizeControlEJ(std::string odesolver, size_t nTimesteps, double timestepSize, double eps,
                    double sc, sgpp::base::ScreenOutput* screen = nullptr, double gamma = 0.5);

  /**
   * Std-Destructor
   */
  virtual ~StepsizeControlEJ();
};

}  // namespace solver
}  // namespace sgpp

#endif /* STEPSIZECONTROLEJ_HPP */
