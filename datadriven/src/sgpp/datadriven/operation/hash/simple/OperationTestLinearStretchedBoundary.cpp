// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/algorithm/test_dataset.hpp>

#include <sgpp/datadriven/operation/hash/simple/OperationTestLinearStretchedBoundary.hpp>


#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace datadriven {

    float_t OperationTestLinearStretchedBoundary::test(SGPP::base::DataVector& alpha, SGPP::base::DataMatrix& data, SGPP::base::DataVector& classes) {
      base::LinearBoundaryBasis<unsigned int, unsigned int> base;
      return test_dataset(this->storage, base, alpha, data, classes);
    }

    float_t OperationTestLinearStretchedBoundary::testMSE(SGPP::base::DataVector& alpha, SGPP::base::DataMatrix& data, SGPP::base::DataVector& refValues) {
      base::LinearBoundaryBasis<unsigned int, unsigned int> base;
      return test_dataset_mse(this->storage, base, alpha, data, refValues);
    }

    float_t OperationTestLinearStretchedBoundary::testWithCharacteristicNumber(SGPP::base::DataVector& alpha, SGPP::base::DataMatrix& data, SGPP::base::DataVector& classes, SGPP::base::DataVector& charaNumbers) {
      base::LinearBoundaryBasis<unsigned int, unsigned int> base;
      return test_datasetWithCharacteristicNumber(this->storage, base, alpha, data, classes, charaNumbers, 0);
    }

    void OperationTestLinearStretchedBoundary::calculateROCcurve(SGPP::base::DataVector& alpha, SGPP::base::DataMatrix& data, SGPP::base::DataVector& classes, SGPP::base::DataVector& thresholds, SGPP::base::DataMatrix& ROC_curve) {
      base::LinearBoundaryBasis<unsigned int, unsigned int> base;
      test_calculateROCcurve(this->storage, base, alpha, data, classes, thresholds, ROC_curve);
    }

  }
}
