/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
//@author Chao qi(qic@in.tum.de)

#ifndef OPERATIONLDLINEAR_HPP
#define OPERATIONLDLINEAR_HPP

#include <sgpp/pde/algorithm/StdUpDown.hpp>

namespace sg {
  namespace finance {

    /**
     *Implements the \f$(x \phi_i(x),\phi_j(x))\f$ operator on linear grids (no boundaries)
     *
     * @version $HEAD$
     */
    class OperationLDLinear: public sg::pde::StdUpDown {
      public:
        /**
         * Constructor
         *
         * @param storage the grid's sg::base::GridStorage object
         */
        OperationLDLinear(sg::base::GridStorage* storage);

        /**
         * Destructor
         */
        virtual ~OperationLDLinear();

      protected:
        /**
         * Up-step in dimension <i>dim</i> for \f$(x \phi_i(x),\phi_j(x))\f$.
         * Applies the up-part of the one-dimensional mass matrix in one dimension.
         * Computes \f[\int_{x=0}^1  x \phi_i(x) \phi_j(x) dx.\f]
         *
         * @param dim dimension in which to apply the up-part
         * @param alpha vector of coefficients
         * @param result vector to store the results in
         */
        virtual void up(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim);

        /**
         * Down-step in dimension <i>dim</i> for \f$(x \phi_i(x),\phi_j(x))\f$.
         * Applies the down-part of the one-dimensional mass matrix in one dimension.
         * Computes \f[\int_{x=0}^1  x \phi_i(x) \phi_j(x) dx.\f]
         *
         * @param dim dimension in which to apply the down-part
         * @param alpha vector of coefficients
         * @param result vector to store the results in
         */
        virtual void down(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim);
    };

  }
}

#endif /* OPERATIONLDLINEAR_HPP */
