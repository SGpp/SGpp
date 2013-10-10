/* ****************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef CONJUGATEGRADIENTSSP_HPP
#define CONJUGATEGRADIENTSSP_HPP

#include "base/operation/OperationMatrixSP.hpp"
#include "base/datatypes/DataVectorSP.hpp"

#include "solver/SLESolverSP.hpp"

#include <iostream>

namespace sg {
  namespace solver {

    class ConjugateGradientsSP : public SLESolverSP {
      public:
        /**
         * Std-Constructor
         */
        ConjugateGradientsSP(size_t imax, float epsilon);

        /**
         * Std-Destructor
         */
        virtual ~ConjugateGradientsSP();

        /**
         * function that defines a solve method for an iterative solver. In contrast to the normal solve routine
         * this method operates on sinlge precision data.
         *
         * @param SystemMatrix reference to an sg::base::OperationMatrix Object that implements the matrix vector multiplication
         * @param alpha the sparse grid's coefficients which have to be determined
         * @param b the right hand side of the system of linear equations
         * @param reuse identifies if the alphas, stored in alpha at calling time, should be reused
         * @param verbose prints information during execution of the solver
         * @param max_threshold additional abort criteria for solver
         */
        void solve(sg::base::OperationMatrixSP& SystemMatrix, sg::base::DataVectorSP& alpha, sg::base::DataVectorSP& b, bool reuse = false, bool verbose = false, float max_threshold = -1.0);

        // Define functions for observer pattern in python

        /**
         * function that signals the start of the CG method (used in python)
         */
        virtual void starting();

        /**
         * function that signals the start of the calculation of the CG method (used in python)
         */
        virtual void calcStarting();

        /**
         * function that signals that one iteration step of the CG method has been completed (used in python)
         */
        virtual void iterationComplete();

        /**
         * function that signals the finish of the cg method (used in python)
         */
        virtual void complete();
    };

  }
}

#endif /* CONJUGATEGRADIENTSSP_HPP */
