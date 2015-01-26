/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Sam Maurus (MA thesis)

#ifndef OPERATIONHESTONZLINEAR_HPP
#define OPERATIONHESTONZLINEAR_HPP

#include "pde/algorithm/UpDownOneOpDim.hpp"

namespace sg {
  namespace finance {

    /**
     * Implements the Heston Z-Operation (corresponds to matrix Z in Master's thesis), that is needed
     * the solve the multidimensional Heston
     * equation, on grids with fix Dirichlet-0-Boundaries.
     *
     * @version $HEAD$
     */
    class OperationHestonZLinear : public sg::pde::UpDownOneOpDim {
      public:
        /**
         * Constructor
         *
         * @param storage the grid's sg::base::GridStorage object
         * @param coef reference to a sg::base::DataVector object that contains the bilinear form's constant coefficients
         */
        OperationHestonZLinear(sg::base::GridStorage* storage, sg::base::DataVector& coef);

        /**
         * Destructor
         */
        virtual ~OperationHestonZLinear();

      protected:

        /**
         * Up-step in dimension <i>dim</i> for \f$(\phi_i(x),\phi_j(x))_{L_2}\f$.
         * Applies the up-part of the one-dimensional mass matrix in one dimension.
         * Computes \f[\int_{x=0}^1  \phi_i(x) \sum_{j, l_i < l_j} \alpha_j \phi_j(x) dx.\f]
         *
         * @param dim dimension in which to apply the up-part
         * @param alpha vector of coefficients
         * @param result vector to store the results in
         */
        virtual void up(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim);

        /**
         * Down-step in dimension <i>dim</i> for \f$(\phi_i(x),\phi_j(x))_{L_2}\f$.
         * Applies the down-part of the one-dimensional mass matrix in one dimension.
         * Computes \f[\int_{x=0}^1  \phi_i(x) \sum_{j, l_i\geq l_j} \alpha_j \phi_j(x) dx.\f]
         *
         * @param dim dimension in which to apply the down-part
         * @param alpha vector of coefficients
         * @param result vector to store the results in
         */
        virtual void down(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim);

        /**
         * down-Gradient step in dimension <i>dim</i> applies the X dphi phi operation
         * in one dimension
         *
         * @param alpha the coefficients of the gridpoints
         * @param result vector with the result of this operation
         * @param dim the dimension in that down-Gradient is applied
         */
        virtual void downOpDim(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim);

        /**
         * up-Gradient step in dimension <i>dim</i> applies the X dphi phi operation
         * in one dimension
         *
         * @param alpha the coefficients of the gridpoints
         * @param result vector with the result of this operation
         * @param dim the dimension in that down-Gradient is applied
         */
        virtual void upOpDim(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim);
    };

  }
}

#endif /* OPERATIONHESTONZLINEAR_HPP */
