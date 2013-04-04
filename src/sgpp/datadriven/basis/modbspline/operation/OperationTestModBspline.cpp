/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Dirk Pflueger (pflueged@in.tum.de)

#include "datadriven/basis/modbspline/operation/OperationTestModBspline.hpp"

#include "base/exception/operation_exception.hpp"
#include "datadriven/algorithm/test_dataset.hpp"

#include "base/datatypes/DataVector.hpp"

namespace sg {
  namespace datadriven {

    double OperationTestModBspline::test(sg::base::DataVector& alpha, sg::base::DataMatrix& data, sg::base::DataVector& classes) {
      return test_dataset(this->storage, base, alpha, data, classes);
    }

    double OperationTestModBspline::testMSE(sg::base::DataVector& alpha, sg::base::DataMatrix& data, sg::base::DataVector& refValues) {
      return test_dataset_mse(this->storage, base, alpha, data, refValues);
    }

    double OperationTestModBspline::testWithCharacteristicNumber(sg::base::DataVector& alpha, sg::base::DataMatrix& data, sg::base::DataVector& classes, sg::base::DataVector& charaNumbers) {
      return test_datasetWithCharacteristicNumber(this->storage, base, alpha, data, classes, charaNumbers, 0.0);
    }

    void OperationTestModBspline::calculateROCcurve(sg::base::DataVector& alpha, sg::base::DataMatrix& data, sg::base::DataVector& classes, sg::base::DataVector& thresholds, sg::base::DataMatrix& ROC_curve) {
      test_calculateROCcurve(this->storage, base, alpha, data, classes, thresholds, ROC_curve);
    }

  }
}
