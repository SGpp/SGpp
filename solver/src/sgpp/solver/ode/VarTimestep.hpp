// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef VARTIMESTEP_HPP
#define VARTIMESTEP_HPP

#include <sgpp/base/application/ScreenOutput.hpp>
#include <sgpp/solver/ODESolver.hpp>
#include <sgpp/solver/ode/StepsizeControl.hpp>
#include <sgpp/globaldef.hpp>

#include <string>

namespace sgpp {
namespace solver {

/**
 * This class implements a step size control using Adams-Bashforth and Crank-Nicolson
 * for solving ordinary partial equations
 *
 */
class VarTimestep : public StepsizeControl {
 protected:
  void predictor(SLESolver& LinearSystemSolver,
                 sgpp::solver::OperationParabolicPDESolverSystem& System, double tmp_timestepsize,
                 sgpp::base::DataVector& dv, sgpp::base::DataVector& corr,
                 sgpp::base::DataVector* rhs);
  void corrector(SLESolver& LinearSystemSolver,
                 sgpp::solver::OperationParabolicPDESolverSystem& System, double tmp_timestepsize,
                 sgpp::base::DataVector& dv, sgpp::base::DataVector* rhs);

  virtual double nextTimestep(double tmp_timestepsize, double tmp_timestepsize_old, double norm,
                               double epsilon);

  std::string _predictor;
  std::string _corrector;

 public:
  /**
   * Std-Constructer
   *
   * @param pred used predictor
   * @param corr used corrector
   * @param imax number of maximum executed iterations
   * @param timestepSize the size of one timestep
   * @param eps the epsilon for the step size control
   * @param screen possible pointer to a sgpp::base::ScreenOutput object
   * @param gamma damping factor
   */
  VarTimestep(std::string pred, std::string corr, size_t imax, double timestepSize, double eps,
              sgpp::base::ScreenOutput* screen = nullptr, double gamma = -1);

  /**
   * Std-Destructor
   */
  virtual ~VarTimestep();
};

}  // namespace solver
}  // namespace sgpp

#endif /* VARTIMESTEP_HPP */
