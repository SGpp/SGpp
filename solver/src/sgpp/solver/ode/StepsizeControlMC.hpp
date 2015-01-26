/* ****************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Peter Hoffmann (peter.hoffmann@mytum.de)

#ifndef STEPSIZECONTROLMC_HPP
#define STEPSIZECONTROLMC_HPP

#include <sgpp/base/application/ScreenOutput.hpp>
#include <sgpp/solver/ode/StepsizeControl.hpp>
#include <string>
#include "VarTimestep.hpp"


namespace sg {
  namespace solver {

    /**
     * This class implements a step size control using Adams-Bashforth and Crank-Nicolson
     * for solving ordinary partial equations
     *
     * @version $HEAD$
     */
    class StepsizeControlMC : public VarTimestep {
      protected:

        double nextTimestep(double tmp_timestepsize, double tmp_timestepsize_old, double norm, double epsilon);

      public:
        /**
         * Std-Constructer
         *
         * @param imax number of maximum executed iterations
         * @param timestepSize the size of one timestep
         * @param eps the epsilon for the step size control
         * @param screen possible pointer to a ScreenOutput object
         */
        StepsizeControlMC(size_t imax, double timestepSize, double eps, sg::base::ScreenOutput* screen = NULL);

        /**
         * Std-Destructor
         */
        virtual ~StepsizeControlMC();
    };

  }
}

#endif /* STEPSIZECONTROLMC_HPP */
