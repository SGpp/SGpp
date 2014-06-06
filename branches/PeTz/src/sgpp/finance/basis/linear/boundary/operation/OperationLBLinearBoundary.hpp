/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Chao qi(qic@in.tum.de)

#ifndef OPERATIONLBLINEARBOUNDARY_HPP
#define OPERATIONLBLINEARBOUNDARY_HPP

#include "pde/algorithm/StdUpDown.hpp"

namespace sg {
  namespace finance {

    /**
     * Implements the standard \f$(\phi_i(x),d\phi_j(x))\f$ operator on linear boundary grids
     *
     * @version $HEAD$
     */
    class OperationLBLinearBoundary: public sg::pde::StdUpDown {
      public:
        /**
         * Constructor
         *
         * @param storage the grid's sg::base::GridStorage object
         */
        OperationLBLinearBoundary(sg::base::GridStorage* storage);

        /**
         * Destructor
         */
        virtual ~OperationLBLinearBoundary();

      protected:
        /**
         * Up-step in dimension <i>dim</i> for \f$(\phi_i(x),d\phi_j(x))\f$.
         * Applies the up-part of the one-dimensional mass matrix in one dimension.
         * Computes \f[\int_{x=0}^1  \phi_i(x)d\phi_j(x) dx.\f]
         *
         * @param dim dimension in which to apply the up-part
         * @param alpha vector of coefficients
         * @param result vector to store the results in
         */
        virtual void up(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim);

        /**
         * Down-step in dimension <i>dim</i> for \f$(\phi_i(x),d\phi_j(x))\f$.
         * Applies the down-part of the one-dimensional mass matrix in one dimension.
         * Computes \f[\int_{x=0}^1  \phi_i(x)d\phi_j(x) dx.\f]
         *
         * @param dim dimension in which to apply the down-part
         * @param alpha vector of coefficients
         * @param result vector to store the results in
         */
        virtual void down(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim);
    };

  }
}

#endif /* OPERATIONLBLINEARBOUNDARY_HPP */
