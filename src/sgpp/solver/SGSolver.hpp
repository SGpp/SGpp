/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef SGSOLVER_HPP
#define SGSOLVER_HPP

#include "solver/TypesSolver.hpp"

#include <cstddef>

namespace sg {
  namespace solver {

    /**
     * Abstract class that defines a solver used in Sparse Grids
     * Applications
     * @todo (pflueged) move public variables to protected scope, add setters and getters
     */
    class SGSolver {
      protected:
        /// Number of Iterations needed for the solve
        size_t nIterations;
        /// Number of maximum iterations for cg
        size_t nMaxIterations;
        /// residuum
        double residuum;
        /// epsilon needed in the, e.g. final error in the iterative solver, or a timestep
        double myEpsilon;

      public:
        /**
         * Std-Constructor
         *
         * @param nMaximumIterations number of maximum executed iterations
         * @param epsilon the final error in the iterative solver, or the size of one timestep
         */
        SGSolver(size_t nMaximumIterations, double epsilon) : nMaxIterations(nMaximumIterations), myEpsilon(epsilon) {
          nIterations = 0;
          residuum = 0.0;
        }

        /**
         * Std-Destructor
         */
        virtual ~SGSolver() { }


        /**
         * function that returns the number of needed solve steps
         *
         * @return the number of needed solve steps of the sovler
         */
        size_t getNumberIterations() {
          return nIterations;
        }

        /**
         * function the returns the residuum (current or final), error of the solver
         *
         * @return the residuum
         */
        double getResiduum() {
          return residuum;
        }

        /**
         * resets the number of maximum iterations
         *
         * @param nIterations the new number of maximum iterations
         */
        void setMaxIterations(size_t nIterations) {
          nMaxIterations = nIterations;
        }

        /**
         * resets the epsilon, that is used in the SGSolver
         *
         * @param eps the new value of epsilon
         */
        void setEpsilon(double eps) {
          myEpsilon = eps;
        }

        /**
         * gets the the epsilon, that is used in the SGSolver
         *
         * @return the epsilon, used in the solver
         */
        double getEpsilon() {
          return myEpsilon;
        }
    };

  }
}

#endif /* SGSOLVER_HPP */
