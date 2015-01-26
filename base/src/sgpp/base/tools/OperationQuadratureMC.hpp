/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Dirk Pflueger (pflueged@in.tum.de)

#ifndef OPERATIONQUADRATUREMC_HPP
#define OPERATIONQUADRATUREMC_HPP

#include <sgpp/base/operation/OperationQuadrature.hpp>
#include <sgpp/base/grid/Grid.hpp>

namespace sg {
  namespace base {

    /**
     * Typedef for general functions that can be passed to integration methods. Requires three parameters. First, the dimensionality, then dim-many coordinates, and then further client data for the function at hand.
     */
    typedef double (*FUNC)(int, double*, void*);


    /**
     * Quadrature on any sparse grid (that has OperationMultipleEval implemented)
     * using Monte Carlo.
     *
     */
    class OperationQuadratureMC : public OperationQuadrature {
      public:
        /**
         * Constructor of OperationQuadratureMC, specifying a grid
         * object and the number of samples to use.
         *
         * @param grid Reference to the grid object
         * @param mcPaths Number of Monte Carlo samples
         *
         * @todo (pflueged) extend to error computation / arbitrary functions
         */
        OperationQuadratureMC(sg::base::Grid& grid, int mcPaths);

        virtual ~OperationQuadratureMC() {}

        /**
         * Quadrature using simple MC in @f$\Omega=[0,1]^d@f$.
         *
         * @param alpha Coefficient vector for current grid
         */
        virtual double doQuadrature(sg::base::DataVector& alpha);

        /**
         * Quadrature of an arbitrary function using
         * simple MC in @f$\Omega=[0,1]^d@f$.
         *
         * @param func The function to integrate
         * @param clientdata Optional data to pass to FUNC
         */
        double doQuadratureFunc(FUNC func, void* clientdata);

        /**
         * Quadrature of the @f$L^2@f$-norm of the error,
         * @f$ ||f(x)-u(x)||_{L^2} @f$, between a given function and the
         * current sparse grid function using
         * simple MC in @f$\Omega=[0,1]^d@f$.
         *
         * @param func The function @f$f(x)@f$
         * @param clientdata Optional data to pass to FUNC
         * @param alpha Coefficient vector for current grid
         */
        double doQuadratureL2Error(FUNC func, void* clientdata, sg::base::DataVector& alpha);

      protected:
        // Pointer to the grid object
        sg::base::Grid* grid;
        // Number of MC paths
        size_t mcPaths;
    };

  }
}

#endif /* OPERATIONQUADRATUREMC_HPP */
