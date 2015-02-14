// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef OPERATIONTESTLINEARSTRETCHED_HPP
#define OPERATIONTESTLINEARSTRETCHED_HPP

#include <sgpp/datadriven/operation/hash/OperationTest.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace datadriven {

    /**
     * This class implements OperationTest for a grids with linearstretched basis ansatzfunctions without boundaries
     */
    class OperationTestLinearStretched : public OperationTest {
      public:
        /**
         * Constructor of OperationTestLinearStretched
         *
         * @param storage Pointer to the grid's gridstorage obejct
         */
        OperationTestLinearStretched(SGPP::base::GridStorage* storage) : storage(storage) {}

        /**
         * Destructor
         */
        virtual ~OperationTestLinearStretched() {}

        virtual float_t test(SGPP::base::DataVector& alpha, SGPP::base::DataMatrix& data, SGPP::base::DataVector& classes);
        virtual float_t testMSE(SGPP::base::DataVector& alpha, SGPP::base::DataMatrix& data, SGPP::base::DataVector& refValues);
        virtual float_t testWithCharacteristicNumber(SGPP::base::DataVector& alpha, SGPP::base::DataMatrix& data, SGPP::base::DataVector& classes, SGPP::base::DataVector& charaNumbers);
        virtual void calculateROCcurve(SGPP::base::DataVector& alpha, SGPP::base::DataMatrix& data, SGPP::base::DataVector& classes, SGPP::base::DataVector& thresholds, SGPP::base::DataMatrix& ROC_curve);

      protected:
        /// Pointer to the grid's GridStorage object
        SGPP::base::GridStorage* storage;
    };

  }
}

#endif /* OPERATIONTESTLINEARSTRETCHED_HPP */