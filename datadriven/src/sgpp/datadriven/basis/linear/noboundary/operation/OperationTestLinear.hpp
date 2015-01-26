/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONTESTLINEAR_HPP
#define OPERATIONTESTLINEAR_HPP

#include <sgpp/datadriven/operation/OperationTest.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
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
        OperationTestLinear(SGPP::base::GridStorage* storage) : storage(storage) {}

        /**
         * Destructor
         */
        virtual ~OperationTestLinear() {}

        virtual double test(SGPP::base::DataVector& alpha, SGPP::base::DataMatrix& data, SGPP::base::DataVector& classes);
        virtual double testMSE(SGPP::base::DataVector& alpha, SGPP::base::DataMatrix& data, SGPP::base::DataVector& refValues);
        virtual double testWithCharacteristicNumber(SGPP::base::DataVector& alpha, SGPP::base::DataMatrix& data, SGPP::base::DataVector& classes, SGPP::base::DataVector& charaNumbers);
        virtual void calculateROCcurve(SGPP::base::DataVector& alpha, SGPP::base::DataMatrix& data, SGPP::base::DataVector& classes, SGPP::base::DataVector& thresholds, SGPP::base::DataMatrix& ROC_curve);

      protected:
        /// Pointer to the grid's GridStorage object
        SGPP::base::GridStorage* storage;
    };

  }
}

#endif /* OPERATIONTESTLINEAR_HPP */
