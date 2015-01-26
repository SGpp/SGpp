// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef OPERATIONLAPLACEENHANCEDLINEARBOUNDARY_HPP
#define OPERATIONLAPLACEENHANCEDLINEARBOUNDARY_HPP

#include <sgpp/misc/pde/algorithm/UpDownOneOpDimEnhanced.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace pde {

    /**
     * Implements the Laplace operator based on
     * the UpDownOneOpDimEnhanced method.
     *
     * @version $HEAD$
     */
    class OperationLaplaceEnhancedLinearBoundary: public UpDownOneOpDimEnhanced {
      public:
        /**
         * Constructor
         *
         * @param storage the grid's SGPP::base::GridStorage object
         */
        OperationLaplaceEnhancedLinearBoundary(SGPP::base::GridStorage* storage);

        /**
         * Constructor of OperationLaplaceLinear
         *
         * @param storage Pointer to the grid's gridstorage obejct
         * @param coef reference to a SGPP::base::DataVector object that contains the bilinear form's constant coefficients; one per dimension
         */
        OperationLaplaceEnhancedLinearBoundary(SGPP::base::GridStorage* storage, SGPP::base::DataVector& coef);

        /**
         * Destructor
         */
        virtual ~OperationLaplaceEnhancedLinearBoundary();

      protected:
        /**
         * Up-step
         *
         * @param dim dimension in which to apply the up-part
         * @param alpha vector of coefficients
         * @param result vector to store the results in
         */
        virtual void up(SGPP::base::DataMatrix& alpha, SGPP::base::DataMatrix& result, size_t dim);

        /**
         * Down-step
         *
         * @param dim dimension in which to apply the down-part
         * @param alpha vector of coefficients
         * @param result vector to store the results in
         */
        virtual void down(SGPP::base::DataMatrix& alpha, SGPP::base::DataMatrix& result, size_t dim);
    };

  }
}

#endif /* OPERATIONLAPLACEENHANCEDLINEARBOUNDARY_HPP */