// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/common/basis/PolyModifiedBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearStretchedBasis.hpp>
#include <sgpp/datadriven/algorithm/test_dataset.hpp>

#include <sgpp/datadriven/operation/hash/simple/OperationTestLinearStretched.hpp>


#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace datadriven {

    float_t OperationTestLinearStretched::test(base::DataVector& alpha, base::DataMatrix& data, base::DataVector& classes) {
      base::LinearStretchedBasis<unsigned int, unsigned int> base;
      return test_dataset(this->storage, base, alpha, data, classes);
    }

    float_t OperationTestLinearStretched::testMSE(base::DataVector& alpha, base::DataMatrix& data, base::DataVector& refValues) {
      base::LinearStretchedBasis<unsigned int, unsigned int> base;
      return test_dataset_mse(this->storage, base, alpha, data, refValues);
    }

    float_t OperationTestLinearStretched::testWithCharacteristicNumber(base::DataVector& alpha, base::DataMatrix& data, base::DataVector& classes, base::DataVector& charaNumbers) {
      base::LinearStretchedBasis<unsigned int, unsigned int> base;
      return test_datasetWithCharacteristicNumber(this->storage, base, alpha, data, classes, charaNumbers, 0.0);
    }

    void OperationTestLinearStretched::calculateROCcurve(SGPP::base::DataVector& alpha, SGPP::base::DataMatrix& data, SGPP::base::DataVector& classes, SGPP::base::DataVector& thresholds, SGPP::base::DataMatrix& ROC_curve) {
      base::LinearStretchedBasis<unsigned int, unsigned int> base;
      test_calculateROCcurve(this->storage, base, alpha, data, classes, thresholds, ROC_curve);
    }

  }
}
