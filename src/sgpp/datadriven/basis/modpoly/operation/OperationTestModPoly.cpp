/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "datadriven/algorithm/test_dataset.hpp"

#include "datadriven/basis/modpoly/operation/OperationTestModPoly.hpp"

#include "base/exception/operation_exception.hpp"

namespace sg {
  namespace datadriven {

    double OperationTestModPoly::test(sg::base::DataVector& alpha, sg::base::DataMatrix& data, sg::base::DataVector& classes) {
      return test_dataset(this->storage, base, alpha, data, classes);
    }

    double OperationTestModPoly::testMSE(sg::base::DataVector& alpha, sg::base::DataMatrix& data, sg::base::DataVector& refValues) {
      return test_dataset_mse(this->storage, base, alpha, data, refValues);
    }

    double OperationTestModPoly::testWithCharacteristicNumber(sg::base::DataVector& alpha, sg::base::DataMatrix& data, sg::base::DataVector& classes, sg::base::DataVector& charaNumbers) {
      return test_datasetWithCharacteristicNumber(this->storage, base, alpha, data, classes, charaNumbers, 0.0);
    }

    void OperationTestModPoly::calculateROCcurve(sg::base::DataVector& alpha, sg::base::DataMatrix& data, sg::base::DataVector& classes, sg::base::DataVector& thresholds, sg::base::DataMatrix& ROC_curve) {
      test_calculateROCcurve(this->storage, base, alpha, data, classes, thresholds, ROC_curve);
    }

  }
}
