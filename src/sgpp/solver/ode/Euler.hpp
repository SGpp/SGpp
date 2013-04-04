/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef EULER_HPP
#define EULER_HPP

#include "base/application/ScreenOutput.hpp"
#include "solver/ODESolver.hpp"
#include <string>

namespace sg {
  namespace solver {

    /**
     * This class implements the explicit and implicit Euler method
     * for solving ordinary partial equations
     *
     * @version $HEAD$
     */
    class Euler : public ODESolver {
      private:
        /// specifies if a grid evaluation should be execute in every time step
        bool bAnimation;
        /// specifies the evaluation per dimension when a animation is created
        size_t evalsAnimation;
        /// specifies the type of euler that should be executed
        std::string ExMode;
        /// Pointer to sg::base::ScreenOutput object
        sg::base::ScreenOutput* myScreen;

      public:
        /**
         * Std-Constructer
         *
         * @param Mode the mode of the euler that should be executed, must be ExEul or ImEul
         * @param imax number of maximum executed iterations
         * @param timestepSize the size of one timestep
         * @param generateAnimation set this, if you want to create a grid evaluation in every time step, in order to create an animation
         * @param numEvalsAnimation specifies the evaluation per dimension when a animation is created
         * @param screen possible pointer to a sg::base::ScreenOutput object
         */
        Euler(std::string Mode, size_t imax, double timestepSize, bool generateAnimation = false, size_t numEvalsAnimation = 20, sg::base::ScreenOutput* screen = NULL);

        /**
         * Std-Destructor
         */
        virtual ~Euler();

        virtual void solve(SLESolver& LinearSystemSolver, sg::pde::OperationParabolicPDESolverSystem& System, bool bIdentifyLastStep = false, bool verbose = false);
    };

  }
}

#endif /* EULER_HPP */
