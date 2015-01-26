/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */

#ifndef STEPSIZECONTROLEJ_HPP
#define STEPSIZECONTROLEJ_HPP

#include <sgpp/base/application/ScreenOutput.hpp>
#include <sgpp/solver/ODESolver.hpp>
#include <sgpp/solver/ode/StepsizeControl.hpp>

namespace sg {
  namespace solver {

    /**
     * This class implements a time step size control using 1D-Diffusion
     * for solving ordinary partial equations
     *
     * For solving the system of linear equations the
     * already implemented CG-method is used
     *
     * @version $HEAD$
     */
    class StepsizeControlEJ : public StepsizeControl {
      private:
        virtual void predictor(SLESolver& LinearSystemSolver, sg::pde::OperationParabolicPDESolverSystem& System,
                               double tmp_timestepsize, sg::base::DataVector& dv, sg::base::DataVector& corr, sg::base::DataVector* rhs);

        virtual void corrector(SLESolver& LinearSystemSolver, sg::pde::OperationParabolicPDESolverSystem& System, double tmp_timestepsize, sg::base::DataVector& dv, sg::base::DataVector* rhs);

        virtual double nextTimestep(double tmp_timestepsize, double tmp_timestepsize_old, double norm, double epsilon);

        virtual double norm(sg::pde::OperationParabolicPDESolverSystem& System, sg::base::DataVector& dv1, sg::base::DataVector& dv2);
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
         * @param screen possible pointer to a sg::base::ScreenOutput object
         * @param gamma used damping factor, default is 0.5
         */
        StepsizeControlEJ(std::string odesolver, size_t nTimesteps, double timestepSize, double eps, double sc, sg::base::ScreenOutput* screen = NULL, double gamma = 0.5);

        /**
         * Std-Destructor
         */
        virtual ~StepsizeControlEJ();
    };

  }
}

#endif /* STEPSIZECONTROLEJ_HPP */
