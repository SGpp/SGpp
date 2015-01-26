/* ****************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Peter Hoffmann (peter.hoffmann@mytum.de)

#ifndef VARTIMESTEP_HPP
#define VARTIMESTEP_HPP

#include <sgpp/base/application/ScreenOutput.hpp>
#include <sgpp/solver/ODESolver.hpp>
#include <sgpp/solver/ode/StepsizeControl.hpp>
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
    class VarTimestep : public StepsizeControl {
      protected:
        void predictor(SLESolver& LinearSystemSolver, sg::pde::OperationParabolicPDESolverSystem& System,
                       double tmp_timestepsize, sg::base::DataVector& dv, sg::base::DataVector& corr, sg::base::DataVector* rhs);
        void corrector(SLESolver& LinearSystemSolver, sg::pde::OperationParabolicPDESolverSystem& System, double tmp_timestepsize, sg::base::DataVector& dv, sg::base::DataVector* rhs);

        virtual double nextTimestep(double tmp_timestepsize, double tmp_timestepsize_old, double norm, double epsilon);

        std::string _predictor;
        std::string _corrector;


      public:
        /**
         * Std-Constructer
         *
         * @param pred used predictor
         * @param corr used corrector
         * @param imax number of maximum executed iterations
         * @param timestepSize the size of one timestep
         * @param eps the epsilon for the step size control
         * @param screen possible pointer to a sg::base::ScreenOutput object
         * @param gamma damping factor
         */
        VarTimestep(std::string pred, std::string corr, size_t imax, double timestepSize, double eps, sg::base::ScreenOutput* screen = NULL, double gamma = -1);

        /**
         * Std-Destructor
         */
        virtual ~VarTimestep();
    };

  }
}

#endif /* VARTIMESTEP_HPP */
