// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONTESTMODWAVELET_HPP
#define OPERATIONTESTMODWAVELET_HPP

#include <sgpp/datadriven/operation/hash/OperationTest.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace datadriven {

    /**
     * This class implements OperationTest for a grid with mod wavelet basis ansatzfunctions
     *
     * @version $HEAD$
     */
    class OperationTestModWavelet : public OperationTest {
      public:
        /**
         * Constructor
         *
         * @param storage the grid's SGPP::base::GridStorage object
         */
        OperationTestModWavelet(SGPP::base::GridStorage* storage) : storage(storage) {}

        /**
         * Destructor
         */
        virtual ~OperationTestModWavelet() {}

        virtual float_t test(SGPP::base::DataVector& alpha, SGPP::base::DataMatrix& data, SGPP::base::DataVector& classes);
        virtual float_t testMSE(SGPP::base::DataVector& alpha, SGPP::base::DataMatrix& data, SGPP::base::DataVector& refValues);
        virtual float_t testWithCharacteristicNumber(SGPP::base::DataVector& alpha, SGPP::base::DataMatrix& data, SGPP::base::DataVector& classes, SGPP::base::DataVector& charaNumbers);
        virtual void calculateROCcurve(SGPP::base::DataVector& alpha, SGPP::base::DataMatrix& data, SGPP::base::DataVector& classes, SGPP::base::DataVector& thresholds, SGPP::base::DataMatrix& ROC_curve);

      protected:
        /// Pointer to SGPP::base::GridStorage object
        SGPP::base::GridStorage* storage;
    };

  }
}

#endif /* OPERATIONTESTMODWAVELET_HPP */