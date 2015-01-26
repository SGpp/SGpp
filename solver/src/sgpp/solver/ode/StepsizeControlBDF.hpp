/* ****************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Peter Hoffmann (peter.hoffmann@mytum.de)


#ifndef STEPSIZECONTROLBDF_HPP
#define STEPSIZECONTROLBDF_HPP

#include <sgpp/base/application/ScreenOutput.hpp>
#include <sgpp/solver/ODESolver.hpp>
#include "VarTimestep.hpp"

namespace sg {
  namespace solver {

    /**
     * This class implements a step size control using the midpoint method and BDF2
     * for solving ordinary partial equations
     *
     * For solving the system of linear equations the
     * already implemented CG-method is used
     *
     * @version $HEAD$
     */
    class StepsizeControlBDF  : public VarTimestep {
      protected:

        double nextTimestep(double tmp_timestepsize, double tmp_timestepsize_old, double norm, double epsilon);

      public:
        /**
         * Std-Constructer
         *
         * @param nTimesteps number of maximum executed iterations
         * @param timestepSize the size of one timestep
         * @param eps the epsilon for the step size control
         * @param screen possible pointer to a sg::base::ScreenOutput object
         */
        StepsizeControlBDF(size_t nTimesteps, double timestepSize, double eps, sg::base::ScreenOutput* screen = NULL);

        /**
         * Std-Destructor
         */
        virtual ~StepsizeControlBDF();

    };

  }
}

#endif /* STEPSIZECONTROLBDF_HPP */
