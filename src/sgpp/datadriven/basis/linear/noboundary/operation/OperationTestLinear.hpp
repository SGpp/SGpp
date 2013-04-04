/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONTESTLINEAR_HPP
#define OPERATIONTESTLINEAR_HPP

#include "datadriven/operation/OperationTest.hpp"
#include "base/grid/GridStorage.hpp"

namespace sg {
  namespace datadriven {

    /**
     * This class implements OperationTest for a grids with linear basis ansatzfunctions without boundaries
     */
    class OperationTestLinear : public OperationTest {
      public:
        /**
         * Constructor of OperationTestLinear
         *
         * @param storage Pointer to the grid's gridstorage obejct
         */
        OperationTestLinear(sg::base::GridStorage* storage) : storage(storage) {}

        /**
         * Destructor
         */
        virtual ~OperationTestLinear() {}

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

#endif /* OPERATIONTESTLINEAR_HPP */
