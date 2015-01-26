/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Dirk Pflueger (pflueged@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONTESTPREWAVELET_HPP
#define OPERATIONTESTPREWAVELET_HPP

#include "datadriven/operation/OperationTest.hpp"
#include "base/grid/GridStorage.hpp"

namespace sg {
  namespace datadriven {

    /**
     * This class implements OperationTest for a grids with prewavelet basis ansatzfunctions without boundaries
     */
    class OperationTestPrewavelet : public OperationTest {
      public:
        /**
         * Constructor of OperationTestPrewavelet
         *
         * @param storage Pointer to the grid's gridstorage obejct
         */
        OperationTestPrewavelet(sg::base::GridStorage* storage) : storage(storage) {}

        /**
         * Destructor
         */
        virtual ~OperationTestPrewavelet() {}

        virtual double test(sg::base::DataVector& alpha, sg::base::DataMatrix& data, sg::base::DataVector& classes);
        virtual double testMSE(sg::base::DataVector& alpha, sg::base::DataMatrix& data, sg::base::DataVector& refValues);
        virtual double testWithCharacteristicNumber(sg::base::DataVector& alpha, sg::base::DataMatrix& data, sg::base::DataVector& classes, sg::base::DataVector& charaNumbers);
        virtual void calculateROCcurve(sg::base::DataVector& alpha, sg::base::DataMatrix& data, sg::base::DataVector& classes, sg::base::DataVector& thresholds, sg::base::DataMatrix& ROC_curve);

      protected:
        /// Pointer to the grid's GridStorage object
        sg::base::GridStorage* storage;
    };

  }
}

#endif /* OPERATIONTESTPREWAVELET_HPP */
