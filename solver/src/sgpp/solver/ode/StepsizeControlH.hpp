/* ****************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Peter Hoffmann (peter.hoffmann@mytum.de)

#ifndef STEPSIZECONTROLH_HPP
#define STEPSIZECONTROLH_HPP

#include <sgpp/base/application/ScreenOutput.hpp>
#include <sgpp/solver/ODESolver.hpp>
#include <string>
#include "StepsizeControl.hpp"

//
#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace solver {

    /**
     * This class implements a step size control using Crank-Nicolson with different step sizes
     * for solving ordinary partial equations
     *
     * @version $HEAD$
     */
    class StepsizeControlH : public StepsizeControl {
      private:


        void predictor(SLESolver& LinearSystemSolver, SGPP::pde::OperationParabolicPDESolverSystem& System,
                       double tmp_timestepsize, SGPP::base::DataVector& dv, SGPP::base::DataVector& corr, SGPP::base::DataVector* rhs);
        void corrector(SLESolver& LinearSystemSolver, SGPP::pde::OperationParabolicPDESolverSystem& System, double tmp_timestepsize, SGPP::base::DataVector& dv, SGPP::base::DataVector* rhs);

        //double twoNorm(SGPP::base::DataVector &dv1, SGPP::base::DataVector &dv2);

        double nextTimestep(double tmp_timestepsize, double tmp_timestepsize_old, double norm, double epsilon);

        std::string _odesolver;


      public:
        /**
         * Std-Constructer
         *
         * @param odesolver the mode of the euler that should be executed, must be ExEul or ImEul
         * @param imax number of maximum executed iterations
         * @param timestepSize the size of one timestep
         * @param eps the epsilon for the step size control
         * @param screen possible pointer to a SGPP::base::ScreenOutput object
         * @param gamma damping factor
         */
        StepsizeControlH(std::string odesolver, size_t imax, double timestepSize, double eps, SGPP::base::ScreenOutput* screen = NULL, double gamma = 0.9);

        /**
         * Std-Destructor
         */
        virtual ~StepsizeControlH();
    };

  }
}

#endif /* STEPSIZECONTROLH_HPP */
