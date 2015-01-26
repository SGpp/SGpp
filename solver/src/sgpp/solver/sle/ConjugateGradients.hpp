/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef CONJUGATEGRADIENTS_HPP
#define CONJUGATEGRADIENTS_HPP

#include <sgpp/solver/SLESolver.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <iostream>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace solver {

    class ConjugateGradients : public SLESolver {
      public:
        /**
         * Std-Constructor
         */
        ConjugateGradients(size_t imax, double epsilon);

        /**
         * Std-Destructor
         */
        virtual ~ConjugateGradients();

        virtual void solve(SGPP::base::OperationMatrix& SystemMatrix, SGPP::base::DataVector& alpha, SGPP::base::DataVector& b, bool reuse = false, bool verbose = false, double max_threshold = -1.0);

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

#endif /* CONJUGATEGRADIENTS_HPP */
