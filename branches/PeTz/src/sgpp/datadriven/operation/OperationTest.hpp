/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONTEST_HPP
#define OPERATIONTEST_HPP

#include "base/datatypes/DataVector.hpp"
#include "base/datatypes/DataMatrix.hpp"

#ifdef _WIN32
#pragma warning(disable: 4267)
#endif

namespace sg {
  namespace datadriven {

    /**
     * Operation the tests the function that is applied the current Sparse sg::base::Grid at a given point
     *
     * @version $HEAD$
     */
    class OperationTest {
      public:
        /**
         * Constructor
         */
        OperationTest() {}

        /**
         * Destructor
         */
        virtual ~OperationTest() {}

        /**
         * Computes the classification accuracy on some test data.
         *
         * The function is evaluated at the given points. Tests on the classes {+1, -1}, cut-off at 0.
         *
         * @param alpha the coefficients of the sparse grid's base functions
         * @param data the coordinates of the evaluation points
         * @param classes sg::base::DataVector holding the class information
         */
        virtual double test(sg::base::DataVector& alpha, sg::base::DataMatrix& data, sg::base::DataVector& classes) = 0;

        /**
         * Computes the regression accuracy on some test data.
         *
         * The function is evaluated at the given points. Calculates the MSE between
         * between the given values and the values evaluated on the sparse grid.
         *
         * @param alpha the coefficients of the sparse grid's base functions
         * @param data the coordinates of the evaluation points
         * @param refValues sg::base::DataVector holding the reference function values
         */
        virtual double testMSE(sg::base::DataVector& alpha, sg::base::DataMatrix& data, sg::base::DataVector& refValues) = 0;

        /**
         * Computes the classification accuracy on some test data.
         *
         * The function is evaluated at the given points. Tests on the classes {+1, -1}, cut-off at 0.
         *
         * Also the number of the TP TN FP FN are determined
         *
         * @param alpha the coefficients of the sparse grid's base functions
         * @param data the coordinates of the evaluation points
         * @param classes sg::base::DataVector the holds the class information
         * @param charaNumbers the number of true positives, true negatives, false positives, false negatives (Vector of length 4)
         */
        virtual double testWithCharacteristicNumber(sg::base::DataVector& alpha, sg::base::DataMatrix& data, sg::base::DataVector& classes, sg::base::DataVector& charaNumbers) = 0;

        /**
         * Computes the classification accuracy on some test data.
         *
         * The function is evaluated at the given points. Tests on the classes {+1, -1}, cut-off variable at threshold
         *
         * Also the number of the TP TN FP FN are determined
         *
         * @param alpha the coefficients of the sparse grid's base functions
         * @param data the coordinates of the evaluation points
         * @param classes sg::base::DataVector the holds the class information
         * @param thresholds the thresholds (between -1.0 and 1.0) for calculating the ROC curve
         * @param ROC_curve DataMatrix into which the ROC curve is stored
         */
        virtual void calculateROCcurve(sg::base::DataVector& alpha, sg::base::DataMatrix& data, sg::base::DataVector& classes, sg::base::DataVector& thresholds, sg::base::DataMatrix& ROC_curve) = 0;
    };

  }
}

#endif /* OPERATIONTEST_HPP */
