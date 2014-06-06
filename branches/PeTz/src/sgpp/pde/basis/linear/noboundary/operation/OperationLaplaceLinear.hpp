/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONLAPLACELINEAR_HPP
#define OPERATIONLAPLACELINEAR_HPP

#include "pde/algorithm/UpDownOneOpDim.hpp"

namespace sg {
  namespace pde {

    /**
     * Implementation for linear functions of Laplace Operation, linear grids without boundaries
     *
     * @version $HEAD$
     */
    class OperationLaplaceLinear: public UpDownOneOpDim {
      public:
        /**
         * Constructor of OperationLaplaceLinear
         *
         * @param storage Pointer to the grid's gridstorage obejct
         */
        OperationLaplaceLinear(sg::base::GridStorage* storage);

        /**
         * Constructor of OperationLaplaceLinear
         *
         * @param storage Pointer to the grid's gridstorage obejct
         * @param coef reference to a sg::base::DataVector object that contains the bilinear form's constant coefficients; one per dimension
         */
        OperationLaplaceLinear(sg::base::GridStorage* storage, sg::base::DataVector& coef);

        /**
         * Destructor
         */
        virtual ~OperationLaplaceLinear();

        virtual void specialOP(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim, size_t gradient_dim);

        virtual void up(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim);

        virtual void down(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim);

        virtual void downOpDim(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim);

        virtual void upOpDim(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim);
    };

  }
}

#endif /* OPERATIONLAPLACELINEAR_HPP */
