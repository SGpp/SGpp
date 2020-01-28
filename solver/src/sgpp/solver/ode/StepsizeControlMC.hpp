// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef STEPSIZECONTROLMC_HPP
#define STEPSIZECONTROLMC_HPP

#include <sgpp/base/application/ScreenOutput.hpp>
#include <sgpp/solver/ode/StepsizeControl.hpp>
#include <sgpp/solver/ode/VarTimestep.hpp>

#include <sgpp/globaldef.hpp>

#include <string>

namespace sgpp {
namespace solver {

/**
 * This class implements a step size control using Adams-Bashforth and Crank-Nicolson
 * for solving ordinary partial equations
 *
 */
class StepsizeControlMC : public VarTimestep {
 protected:
  double nextTimestep(double tmp_timestepsize, double tmp_timestepsize_old, double norm,
                       double epsilon);

 public:
  /**
   * Std-Constructer
   *
   * @param imax number of maximum executed iterations
   * @param timestepSize the size of one timestep
   * @param eps the epsilon for the step size control
   * @param screen possible pointer to a ScreenOutput object
   */
  StepsizeControlMC(size_t imax, double timestepSize, double eps,
                    sgpp::base::ScreenOutput* screen = nullptr);

  /**
   * Std-Destructor
   */
  virtual ~StepsizeControlMC();
};

}  // namespace solver
}  // namespace sgpp

#endif /* STEPSIZECONTROLMC_HPP */
