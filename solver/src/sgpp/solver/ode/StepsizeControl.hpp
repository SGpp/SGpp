// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef STEPSIZECONTROL_HPP
#define STEPSIZECONTROL_HPP

#include <sgpp/base/application/ScreenOutput.hpp>
#include <sgpp/solver/ODESolver.hpp>
#include <sgpp/pde/operation/hash/OperationParabolicPDESolverSystem.hpp>
#include <string>
//
#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace solver {

    /**
     * This class implements a step size control using Adams-Bashforth and Crank-Nicolson
     * for solving ordinary partial equations
     *
     * @version $HEAD$
     */
    class StepsizeControl : public ODESolver {
      private:
        /// identify of grid coarsening is used
        bool useCoarsen;


      protected:
        /// Pointer to SGPP::base::ScreenOutput object
        SGPP::base::ScreenOutput* myScreen;

        /// temp. Stepsize Control
        float_t mySC;

        /// epsilon for the step size control
        float_t myEps;

        virtual void predictor(SLESolver& LinearSystemSolver, SGPP::pde::OperationParabolicPDESolverSystem& System,
                               float_t tmp_timestepsize, SGPP::base::DataVector& dv, SGPP::base::DataVector& corr, SGPP::base::DataVector* rhs) = 0;
        virtual void corrector(SLESolver& LinearSystemSolver, SGPP::pde::OperationParabolicPDESolverSystem& System, float_t tmp_timestepsize, SGPP::base::DataVector& dv, SGPP::base::DataVector* rhs) = 0;

        virtual float_t norm(SGPP::pde::OperationParabolicPDESolverSystem& System, SGPP::base::DataVector& dv1, SGPP::base::DataVector& dv2);

        virtual float_t nextTimestep(float_t tmp_timestepsize, float_t tmp_timestepsize_old, float_t norm, float_t epsilon) = 0;

        float_t twoNorm(SGPP::pde::OperationParabolicPDESolverSystem& System, SGPP::base::DataVector& dv1, SGPP::base::DataVector& dv2);

        float_t maxNorm(SGPP::pde::OperationParabolicPDESolverSystem& System, SGPP::base::DataVector& dv1, SGPP::base::DataVector& dv2);

        std::string filename;

        /// damping factor
        float_t _gamma;


      public:
        /**
         * Std-Constructer
         *
         * @param sc step size
         * @param imax number of maximum executed iterations
         * @param timestepSize the size of one timestep
         * @param eps the epsilon for the step size control
         * @param screen possible pointer to a SGPP::base::ScreenOutput object
         * @param gamma damping factor
         */
        StepsizeControl(size_t imax, float_t timestepSize, float_t eps, float_t sc, SGPP::base::ScreenOutput* screen = NULL, float_t gamma = 0.5);

        /**
         * Std-Destructor
         */
        virtual ~StepsizeControl();

        void solve(SLESolver& LinearSystemSolver, SGPP::pde::OperationParabolicPDESolverSystem& System, bool bIdentifyLastStep = false, bool verbose = false);
    };

  }
}

#endif /* STEPSIZECONTROL_H_ */