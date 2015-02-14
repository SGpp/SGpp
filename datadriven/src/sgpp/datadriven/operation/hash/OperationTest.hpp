// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef OPERATIONTEST_HPP
#define OPERATIONTEST_HPP

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>

#ifdef _WIN32
#pragma warning(disable: 4267)
#endif

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace datadriven {

    /**
     * Operation the tests the function that is applied the current Sparse SGPP::base::Grid at a given point
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
         * @param classes SGPP::base::DataVector holding the class information
         */
        virtual float_t test(SGPP::base::DataVector& alpha, SGPP::base::DataMatrix& data, SGPP::base::DataVector& classes) = 0;

        /**
         * Computes the regression accuracy on some test data.
         *
         * The function is evaluated at the given points. Calculates the MSE between
         * between the given values and the values evaluated on the sparse grid.
         *
         * @param alpha the coefficients of the sparse grid's base functions
         * @param data the coordinates of the evaluation points
         * @param refValues SGPP::base::DataVector holding the reference function values
         */
        virtual float_t testMSE(SGPP::base::DataVector& alpha, SGPP::base::DataMatrix& data, SGPP::base::DataVector& refValues) = 0;

        /**
         * Computes the classification accuracy on some test data.
         *
         * The function is evaluated at the given points. Tests on the classes {+1, -1}, cut-off at 0.
         *
         * Also the number of the TP TN FP FN are determined
         *
         * @param alpha the coefficients of the sparse grid's base functions
         * @param data the coordinates of the evaluation points
         * @param classes SGPP::base::DataVector the holds the class information
         * @param charaNumbers the number of true positives, true negatives, false positives, false negatives (Vector of length 4)
         */
        virtual float_t testWithCharacteristicNumber(SGPP::base::DataVector& alpha, SGPP::base::DataMatrix& data, SGPP::base::DataVector& classes, SGPP::base::DataVector& charaNumbers) = 0;

        /**
         * Computes the classification accuracy on some test data.
         *
         * The function is evaluated at the given points. Tests on the classes {+1, -1}, cut-off variable at threshold
         *
         * Also the number of the TP TN FP FN are determined
         *
         * @param alpha the coefficients of the sparse grid's base functions
         * @param data the coordinates of the evaluation points
         * @param classes SGPP::base::DataVector the holds the class information
         * @param thresholds the thresholds (between -1.0 and 1.0) for calculating the ROC curve
         * @param ROC_curve DataMatrix into which the ROC curve is stored
         */
        virtual void calculateROCcurve(SGPP::base::DataVector& alpha, SGPP::base::DataMatrix& data, SGPP::base::DataVector& classes, SGPP::base::DataVector& thresholds, SGPP::base::DataMatrix& ROC_curve) = 0;
    };

  }
}

#endif /* OPERATIONTEST_HPP */