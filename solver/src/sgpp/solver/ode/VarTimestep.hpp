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


namespace SGPP {
namespace solver {

/**
 * This class implements a step size control using Adams-Bashforth and Crank-Nicolson
 * for solving ordinary partial equations
 *
 */
class VarTimestep : public StepsizeControl {
 protected:
  void predictor(SLESolver& LinearSystemSolver,
                 SGPP::solver::OperationParabolicPDESolverSystem& System,
                 float_t tmp_timestepsize, SGPP::base::DataVector& dv,
                 SGPP::base::DataVector& corr, SGPP::base::DataVector* rhs);
  void corrector(SLESolver& LinearSystemSolver,
                 SGPP::solver::OperationParabolicPDESolverSystem& System,
                 float_t tmp_timestepsize, SGPP::base::DataVector& dv,
                 SGPP::base::DataVector* rhs);

  virtual float_t nextTimestep(float_t tmp_timestepsize,
                               float_t tmp_timestepsize_old, float_t norm, float_t epsilon);

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
   * @param screen possible pointer to a SGPP::base::ScreenOutput object
   * @param gamma damping factor
   */
  VarTimestep(std::string pred, std::string corr, size_t imax,
              float_t timestepSize, float_t eps, SGPP::base::ScreenOutput* screen = NULL,
              float_t gamma = -1);

  /**
   * Std-Destructor
   */
  virtual ~VarTimestep();
};

}  // namespace solver
}  // namespace SGPP

#endif /* VARTIMESTEP_HPP */
