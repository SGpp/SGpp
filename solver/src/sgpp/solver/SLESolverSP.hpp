/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef SLESOLVERSP_HPP
#define SLESOLVERSP_HPP

#include <sgpp/base/datatypes/DataVectorSP.hpp>
#include <sgpp/base/operation/OperationMatrixSP.hpp>

#include <sgpp/solver/SGSolverSP.hpp>

#ifndef DEFAULT_RES_THRESHOLD
#define DEFAULT_RES_THRESHOLD -1.0
#endif

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace solver {

    class SLESolverSP : public SGSolverSP {
      public:
        /**
         * Std-Constructor
         *
         * @param imax number of maximum executed iterations
         * @param epsilon the final error in the iterative solver
         */
        SLESolverSP(size_t imax, float epsilon) : SGSolverSP(imax, epsilon) {
        }

        /**
         * Std-Destructor
         */
        virtual ~SLESolverSP() { }

        /**
         * Pure virtual Function that defines a solve method for an iterative solver
         *
         * @param SystemMatrix reference to an SGPP::base::OperationMatrix Object that implements the matrix vector multiplication
         * @param alpha the sparse grid's coefficients which have to be determined
         * @param b the right hand side of the system of linear equations
         * @param reuse identifies if the alphas, stored in alpha at calling time, should be reused
         * @param verbose prints information during execution of the solver
         * @param max_threshold additional abort criteria for solver, default value is 10^-9!
         */
        virtual void solve(SGPP::base::OperationMatrixSP& SystemMatrix, SGPP::base::DataVectorSP& alpha, SGPP::base::DataVectorSP& b, bool reuse = false, bool verbose = false, float max_threshold = DEFAULT_RES_THRESHOLD) = 0;
    };

  }
}

#endif /* SLESOLVER_HPP */
