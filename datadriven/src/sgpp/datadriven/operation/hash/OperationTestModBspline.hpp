// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef OPERATIONTESTMODBSPLINE_HPP
#define OPERATIONTESTMODBSPLINE_HPP

#include <sgpp/datadriven/operation/hash/OperationTest.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

#include <sgpp/base/basis/modbspline/ModifiedBsplineBasis.hpp>


#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace datadriven {

    /**
     * This class implements OperationTest for a grids with modified bspline basis functions with a certain degree
     *
     * @version $HEAD$
     */
    class OperationTestModBspline : public OperationTest {
      public:
        /**
         * Constructor
         *
         * @param storage the grid's base::GridStorage object
         * @param degree the bspline's degree
         */
        OperationTestModBspline(base::GridStorage* storage, size_t degree) : storage(storage), base(degree) {}

        /**
         * Destructor
         */
        virtual ~OperationTestModBspline() {}

        virtual double test(base::DataVector& alpha, base::DataMatrix& data, base::DataVector& classes);
        virtual double testMSE(base::DataVector& alpha, base::DataMatrix& data, base::DataVector& refValues);
        virtual double testWithCharacteristicNumber(SGPP::base::DataVector& alpha, SGPP::base::DataMatrix& data, SGPP::base::DataVector& classes, SGPP::base::DataVector& charaNumbers);
        virtual void calculateROCcurve(SGPP::base::DataVector& alpha, SGPP::base::DataMatrix& data, SGPP::base::DataVector& classes, SGPP::base::DataVector& thresholds, SGPP::base::DataMatrix& ROC_curve);

      protected:
        /// Pointer to base::GridStorage object
        base::GridStorage* storage;
        /// Mod Bspline Basis object
        base::SModBsplineBase base;
    };

  }
}

#endif /* OPERATIONTESTMODBSPLINE_HPP */