// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/operation/hash/simple/OperationTestModBspline.hpp>

#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/datadriven/algorithm/test_dataset.hpp>

#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {

double OperationTestModBspline::test(sgpp::base::DataVector& alpha, sgpp::base::DataMatrix& data,
                                      sgpp::base::DataVector& classes) {
  return test_dataset(this->storage, base, alpha, data, classes);
}

double OperationTestModBspline::testMSE(sgpp::base::DataVector& alpha,
                                         sgpp::base::DataMatrix& data,
                                         sgpp::base::DataVector& refValues) {
  return test_dataset_mse(this->storage, base, alpha, data, refValues);
}

double OperationTestModBspline::testWithCharacteristicNumber(
    sgpp::base::DataVector& alpha, sgpp::base::DataMatrix& data, sgpp::base::DataVector& classes,
    sgpp::base::DataVector& charaNumbers) {
  return test_datasetWithCharacteristicNumber(this->storage, base, alpha, data, classes,
                                              charaNumbers, 0.0);
}

void OperationTestModBspline::calculateROCcurve(sgpp::base::DataVector& alpha,
                                                sgpp::base::DataMatrix& data,
                                                sgpp::base::DataVector& classes,
                                                sgpp::base::DataVector& thresholds,
                                                sgpp::base::DataMatrix& ROC_curve) {
  test_calculateROCcurve(this->storage, base, alpha, data, classes, thresholds, ROC_curve);
}

}  // namespace datadriven
}  // namespace sgpp
