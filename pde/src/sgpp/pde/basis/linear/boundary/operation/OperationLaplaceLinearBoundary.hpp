/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONLAPLACELINEARBOUNDARY_HPP
#define OPERATIONLAPLACELINEARBOUNDARY_HPP

#include <sgpp/pde/algorithm/UpDownOneOpDim.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace pde {

    /**
     * Implementation of Laplace for linear functions with boundaries
     *
     * @version $HEAD$
     */
    class OperationLaplaceLinearBoundary: public UpDownOneOpDim {
      public:
        /**
         * Constructor
         *
         * @param storage the grid's SGPP::base::GridStorage object
         */
        OperationLaplaceLinearBoundary(SGPP::base::GridStorage* storage);

        /**
         * Constructor of OperationLaplaceLinear
         *
         * @param storage Pointer to the grid's gridstorage obejct
         * @param coef reference to a SGPP::base::DataVector object that contains the bilinear form's constant coefficients; one per dimension
         */
        OperationLaplaceLinearBoundary(SGPP::base::GridStorage* storage, SGPP::base::DataVector& coef);

        /**
         * Destructor
         */
        virtual ~OperationLaplaceLinearBoundary();


      protected:
        virtual void up(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim);

        virtual void down(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim);

        virtual void downOpDim(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim);

        virtual void upOpDim(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim);
    };

  }
}

#endif /* OPERATIONLAPLACELINEARBOUNDARY_HPP */
