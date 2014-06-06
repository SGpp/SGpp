/* ****************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Peter Hoffmann (peter.hoffmann@mytum.de)

#ifndef STEPSIZECONTROL_HPP
#define STEPSIZECONTROL_HPP

#include "base/application/ScreenOutput.hpp"
#include "solver/ODESolver.hpp"
#include "pde/operation/OperationParabolicPDESolverSystem.hpp"
#include <string>
//
namespace sg {
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
        /// Pointer to sg::base::ScreenOutput object
        sg::base::ScreenOutput* myScreen;

        /// temp. Stepsize Control
        double mySC;

        /// epsilon for the step size control
        double myEps;

        virtual void predictor(SLESolver& LinearSystemSolver, sg::pde::OperationParabolicPDESolverSystem& System,
                               double tmp_timestepsize, sg::base::DataVector& dv, sg::base::DataVector& corr, sg::base::DataVector* rhs) = 0;
        virtual void corrector(SLESolver& LinearSystemSolver, sg::pde::OperationParabolicPDESolverSystem& System, double tmp_timestepsize, sg::base::DataVector& dv, sg::base::DataVector* rhs) = 0;

        virtual double norm(sg::pde::OperationParabolicPDESolverSystem& System, sg::base::DataVector& dv1, sg::base::DataVector& dv2);

        virtual double nextTimestep(double tmp_timestepsize, double tmp_timestepsize_old, double norm, double epsilon) = 0;

        double twoNorm(sg::pde::OperationParabolicPDESolverSystem& System, sg::base::DataVector& dv1, sg::base::DataVector& dv2);

        double maxNorm(sg::pde::OperationParabolicPDESolverSystem& System, sg::base::DataVector& dv1, sg::base::DataVector& dv2);

        std::string filename;

        /// damping factor
        double _gamma;


      public:
        /**
         * Std-Constructer
         *
         * @param sc step size
         * @param imax number of maximum executed iterations
         * @param timestepSize the size of one timestep
         * @param eps the epsilon for the step size control
         * @param screen possible pointer to a sg::base::ScreenOutput object
         * @param gamma damping factor
         */
        StepsizeControl(size_t imax, double timestepSize, double eps, double sc, sg::base::ScreenOutput* screen = NULL, double gamma = 0.5);

        /**
         * Std-Destructor
         */
        virtual ~StepsizeControl();

        void solve(SLESolver& LinearSystemSolver, sg::pde::OperationParabolicPDESolverSystem& System, bool bIdentifyLastStep = false, bool verbose = false);
    };

  }
}

#endif /* STEPSIZECONTROL_H_ */
