/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef ODESOLVER_HPP
#define ODESOLVER_HPP

#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/pde/operation/OperationParabolicPDESolverSystem.hpp>

#include <sgpp/solver/SGSolver.hpp>
#include <sgpp/solver/SLESolver.hpp>

namespace sg {
  namespace solver {

    class ODESolver : public SGSolver {
      public:
        /**
         * Std-Constructor
         *
         * @param imax number of maximum executed iterations
         * @param timestepSize the size of one timestep
         */
        ODESolver(size_t imax, double timestepSize) : SGSolver(imax, timestepSize) {
        }

        /**
         * Std-Destructor
         */
        virtual ~ODESolver() { }

        /**
         * Pure virtual Function that defines a solve method for an ODE solver
         *
         * @param LinearSystemSolver reference to an instance of a linear system solver that is used by this ODE solver
         * @param System reference to an sg::base::OperationMatrix Object that implements the matrix vector multiplication
         * @param bIdentifyLastStep set this to true to tell System the last step
         * @param verbose prints information during execution of the solver
         */
        virtual void solve(SLESolver& LinearSystemSolver, sg::pde::OperationParabolicPDESolverSystem& System, bool bIdentifyLastStep = false, bool verbose = false) = 0;
    };

  }
}

#endif /* ODESOLVER_HPP */
