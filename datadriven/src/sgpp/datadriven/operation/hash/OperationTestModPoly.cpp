// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#include <sgpp/datadriven/algorithm/test_dataset.hpp>

#include <sgpp/datadriven/operation/hash/OperationTestModPoly.hpp>

#include <sgpp/base/exception/operation_exception.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace datadriven {

    double OperationTestModPoly::test(SGPP::base::DataVector& alpha, SGPP::base::DataMatrix& data, SGPP::base::DataVector& classes) {
      return test_dataset(this->storage, base, alpha, data, classes);
    }

    double OperationTestModPoly::testMSE(SGPP::base::DataVector& alpha, SGPP::base::DataMatrix& data, SGPP::base::DataVector& refValues) {
      return test_dataset_mse(this->storage, base, alpha, data, refValues);
    }

    double OperationTestModPoly::testWithCharacteristicNumber(SGPP::base::DataVector& alpha, SGPP::base::DataMatrix& data, SGPP::base::DataVector& classes, SGPP::base::DataVector& charaNumbers) {
      return test_datasetWithCharacteristicNumber(this->storage, base, alpha, data, classes, charaNumbers, 0.0);
    }

    void OperationTestModPoly::calculateROCcurve(SGPP::base::DataVector& alpha, SGPP::base::DataMatrix& data, SGPP::base::DataVector& classes, SGPP::base::DataVector& thresholds, SGPP::base::DataMatrix& ROC_curve) {
      test_calculateROCcurve(this->storage, base, alpha, data, classes, thresholds, ROC_curve);
    }

  }
}