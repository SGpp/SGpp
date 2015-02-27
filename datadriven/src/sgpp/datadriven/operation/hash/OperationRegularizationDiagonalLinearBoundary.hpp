// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONREGULARIZATIONDIAGONALLINEARBOUNDARY_HPP
#define OPERATIONREGULARIZATIONDIAGONALLINEARBOUNDARY_HPP

#include <sgpp/datadriven/operation/hash/OperationRegularizationDiagonal.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace datadriven {

    /**
     * Implementation of the application of a diagonal matrix to a
     * DataVector for regularization.
     * This class implements several scaling possibilities for piecewise
     * linear basis functions with and without boundaries.
     */
    class OperationRegularizationDiagonalLinearBoundary: public OperationRegularizationDiagonal {
      protected:
        /**
         * Initialize Hkmix
         * @param k Parameter k
         */
        virtual void initHkmix(float_t k);

        /**
         * Initialize H0HkLaplace
         * @param k Parameter k
         */
        virtual void initH0HkLaplace(float_t k);

      public:

        /**
         * Constructor of OperationRegularizationDiagonalLinearBoundary.
         * @param storage Pointer to grid's storage object
         * @param mode Mode, specifying which regularization to use. Example: OperationRegularizationDiagonal::HKMIX.
         * @param k Parameter for @f$H^k@f$
         */
        OperationRegularizationDiagonalLinearBoundary(base::GridStorage* storage, int mode, float_t k);
        //      : OperationRegularizationDiagonal(storage, mode, k) {
        //      init();
        //    };

    };

  }
}
#endif /* OPERATIONREGULARIZATIONDIAGONALLINEARBOUNDARY_HPP */