/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONTESTMODPOLY_HPP
#define OPERATIONTESTMODPOLY_HPP

#include "datadriven/operation/OperationTest.hpp"
#include "base/grid/GridStorage.hpp"

#include "base/basis/modpoly/ModifiedPolyBasis.hpp"

#include "base/datatypes/DataVector.hpp"

namespace sg {
  namespace datadriven {

    /**
     * This class implements OperationTest for a grids with mod poly basis ansatzfunctions with
     *
     * @version $HEAD$
     */
    class OperationTestModPoly : public OperationTest {
      public:
        /**
         * Constructor
         *
         * @param storage the grid's base::GridStorage object
         * @param degree the polynom's max. degree
         */
        OperationTestModPoly(base::GridStorage* storage, size_t degree) : storage(storage), base(degree) {}

        /**
         * Destructor
         */
        virtual ~OperationTestModPoly() {}

        virtual double test(base::DataVector& alpha, base::DataMatrix& data, base::DataVector& classes);
        virtual double testMSE(base::DataVector& alpha, base::DataMatrix& data, base::DataVector& refValues);
        virtual double testWithCharacteristicNumber(sg::base::DataVector& alpha, sg::base::DataMatrix& data, sg::base::DataVector& classes, sg::base::DataVector& charaNumbers);
        virtual void calculateROCcurve(sg::base::DataVector& alpha, sg::base::DataMatrix& data, sg::base::DataVector& classes, sg::base::DataVector& thresholds, sg::base::DataMatrix& ROC_curve);

      protected:
        /// Pointer to base::GridStorage object
        base::GridStorage* storage;
        /// Mod Poly Basis object
        base::SModPolyBase base;
    };

  }
}

#endif /* OPERATIONTESTMODPOLY_HPP */
