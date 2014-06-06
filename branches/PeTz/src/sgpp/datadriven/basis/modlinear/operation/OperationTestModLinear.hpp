/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONTESTMODLINEAR_HPP
#define OPERATIONTESTMODLINEAR_HPP

#include "datadriven/operation/OperationTest.hpp"
#include "base/grid/GridStorage.hpp"

namespace sg {
  namespace datadriven {

    /**
     * This class implements sg::base::OperationEval for a grids with mod linear basis ansatzfunctions with
     *
     * @version $HEAD$
     */
    class OperationTestModLinear : public OperationTest {
      public:
        /**
         * Constructor
         *
         * @param storage the grid's sg::base::GridStorage object
         */
        OperationTestModLinear(sg::base::GridStorage* storage) : storage(storage) {}

        /**
         * Destructor
         */
        virtual ~OperationTestModLinear() {}

        virtual double test(sg::base::DataVector& alpha, sg::base::DataMatrix& data, sg::base::DataVector& classes);
        virtual double testMSE(sg::base::DataVector& alpha, sg::base::DataMatrix& data, sg::base::DataVector& refValues);
        virtual double testWithCharacteristicNumber(sg::base::DataVector& alpha, sg::base::DataMatrix& data, sg::base::DataVector& classes, sg::base::DataVector& charaNumbers);
        virtual void calculateROCcurve(sg::base::DataVector& alpha, sg::base::DataMatrix& data, sg::base::DataVector& classes, sg::base::DataVector& thresholds, sg::base::DataMatrix& ROC_curve);

      protected:
        /// Pointer to sg::base::GridStorage object
        sg::base::GridStorage* storage;

    };

  }
}

#endif /* OPERATIONTESTMODLINEAR_HPP */
