/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONLAPLACELINEAR_HPP
#define OPERATIONLAPLACELINEAR_HPP

#include <sgpp/pde/algorithm/UpDownOneOpDim.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
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
        OperationLaplaceLinear(SGPP::base::GridStorage* storage);

        /**
         * Constructor of OperationLaplaceLinear
         *
         * @param storage Pointer to the grid's gridstorage obejct
         * @param coef reference to a SGPP::base::DataVector object that contains the bilinear form's constant coefficients; one per dimension
         */
        OperationLaplaceLinear(SGPP::base::GridStorage* storage, SGPP::base::DataVector& coef);

        /**
         * Destructor
         */
        virtual ~OperationLaplaceLinear();

        virtual void specialOP(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim, size_t gradient_dim);

        virtual void up(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim);

        virtual void down(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim);

        virtual void downOpDim(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim);

        virtual void upOpDim(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim);
    };

  }
}

#endif /* OPERATIONLAPLACELINEAR_HPP */
