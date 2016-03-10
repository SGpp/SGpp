// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/operation/hash/simple/OperationTestPoly.hpp>

#include <sgpp/datadriven/algorithm/test_dataset.hpp>

#include <sgpp/base/exception/operation_exception.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {

double OperationTestPoly::test(base::DataVector& alpha, base::DataMatrix& data,
                                base::DataVector& classes) {
  return test_dataset(this->storage, base, alpha, data, classes);
}

double OperationTestPoly::testMSE(base::DataVector& alpha, base::DataMatrix& data,
                                   base::DataVector& refValues) {
  return test_dataset_mse(this->storage, base, alpha, data, refValues);
}

double OperationTestPoly::testWithCharacteristicNumber(base::DataVector& alpha,
                                                        base::DataMatrix& data,
                                                        base::DataVector& classes,
                                                        base::DataVector& charaNumbers) {
  return test_datasetWithCharacteristicNumber(this->storage, base, alpha, data, classes,
                                              charaNumbers, 0.0);
}

void OperationTestPoly::calculateROCcurve(sgpp::base::DataVector& alpha,
                                          sgpp::base::DataMatrix& data,
                                          sgpp::base::DataVector& classes,
                                          sgpp::base::DataVector& thresholds,
                                          sgpp::base::DataMatrix& ROC_curve) {
  test_calculateROCcurve(this->storage, base, alpha, data, classes, thresholds, ROC_curve);
}

}  // namespace datadriven
}  // namespace sgpp
