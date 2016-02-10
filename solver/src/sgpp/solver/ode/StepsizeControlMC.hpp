// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef STEPSIZECONTROLMC_HPP
#define STEPSIZECONTROLMC_HPP

#include <sgpp/base/application/ScreenOutput.hpp>
#include <sgpp/solver/ode/StepsizeControl.hpp>
#include <string>
#include <sgpp/solver/ode/VarTimestep.hpp>


#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace solver {

/**
 * This class implements a step size control using Adams-Bashforth and Crank-Nicolson
 * for solving ordinary partial equations
 *
 */
class StepsizeControlMC : public VarTimestep {
 protected:

  float_t nextTimestep(float_t tmp_timestepsize, float_t tmp_timestepsize_old,
                       float_t norm, float_t epsilon);

 public:
  /**
   * Std-Constructer
   *
   * @param imax number of maximum executed iterations
   * @param timestepSize the size of one timestep
   * @param eps the epsilon for the step size control
   * @param screen possible pointer to a ScreenOutput object
   */
  StepsizeControlMC(size_t imax, float_t timestepSize, float_t eps,
                    SGPP::base::ScreenOutput* screen = NULL);

  /**
   * Std-Destructor
   */
  virtual ~StepsizeControlMC();
};

}
}

#endif /* STEPSIZECONTROLMC_HPP */