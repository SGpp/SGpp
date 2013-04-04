/* ****************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Peter Hoffmann (peter.hoffmann@mytum.de)

#ifndef ADAMSBASHFORTH_HPP
#define ADAMSBASHFORTH_HPP

#include "base/application/ScreenOutput.hpp"
#include "solver/ODESolver.hpp"
#include <string>

namespace sg {
  namespace solver {

    /**
     * This class implements the explicit Adams-Bashforth method
     * for solving ordinary partial equations
     *
     * @version $HEAD$
     */
    class AdamsBashforth : public ODESolver {
      private:
        /// Pointer to sg::base::ScreenOutput object
        sg::base::ScreenOutput* myScreen;

      public:
        /**
         * Std-Constructer
         *
         * @param imax number of maximum executed iterations
         * @param timestepSize the size of one timestep
         * @param screen possible pointer to a sg::base::ScreenOutput object
         */
        AdamsBashforth(size_t imax, double timestepSize, sg::base::ScreenOutput* screen = NULL);

        /**
         * Std-Destructor
         */
        virtual ~AdamsBashforth();

        virtual void solve(SLESolver& LinearSystemSolver, sg::pde::OperationParabolicPDESolverSystem& System, bool bIdentifyLastStep = false, bool verbose = false);
    };

  }
}

#endif /* ADAMSBASHFORTH_HPP */
