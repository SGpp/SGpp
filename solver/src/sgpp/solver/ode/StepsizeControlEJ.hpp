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


namespace SGPP {
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
        virtual void predictor(SLESolver& LinearSystemSolver, SGPP::pde::OperationParabolicPDESolverSystem& System,
                               float_t tmp_timestepsize, SGPP::base::DataVector& dv, SGPP::base::DataVector& corr, SGPP::base::DataVector* rhs);

        virtual void corrector(SLESolver& LinearSystemSolver, SGPP::pde::OperationParabolicPDESolverSystem& System, float_t tmp_timestepsize, SGPP::base::DataVector& dv, SGPP::base::DataVector* rhs);

        virtual float_t nextTimestep(float_t tmp_timestepsize, float_t tmp_timestepsize_old, float_t norm, float_t epsilon);

        virtual float_t norm(SGPP::pde::OperationParabolicPDESolverSystem& System, SGPP::base::DataVector& dv1, SGPP::base::DataVector& dv2);
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
         * @param screen possible pointer to a SGPP::base::ScreenOutput object
         * @param gamma used damping factor, default is 0.5
         */
        StepsizeControlEJ(std::string odesolver, size_t nTimesteps, float_t timestepSize, float_t eps, float_t sc, SGPP::base::ScreenOutput* screen = NULL, float_t gamma = 0.5);

        /**
         * Std-Destructor
         */
        virtual ~StepsizeControlEJ();
    };

  }
}

#endif /* STEPSIZECONTROLEJ_HPP */
