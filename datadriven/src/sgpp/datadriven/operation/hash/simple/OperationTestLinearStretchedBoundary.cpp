// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/algorithm/test_dataset.hpp>

#include <sgpp/datadriven/operation/hash/simple/OperationTestLinearStretchedBoundary.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {

double OperationTestLinearStretchedBoundary::test(sgpp::base::DataVector& alpha,
                                                   sgpp::base::DataMatrix& data,
                                                   sgpp::base::DataVector& classes) {
  base::LinearBoundaryBasis<unsigned int, unsigned int> base;
  return test_dataset(this->storage, base, alpha, data, classes);
}

double OperationTestLinearStretchedBoundary::testMSE(sgpp::base::DataVector& alpha,
                                                      sgpp::base::DataMatrix& data,
                                                      sgpp::base::DataVector& refValues) {
  base::LinearBoundaryBasis<unsigned int, unsigned int> base;
  return test_dataset_mse(this->storage, base, alpha, data, refValues);
}

double OperationTestLinearStretchedBoundary::testWithCharacteristicNumber(
    sgpp::base::DataVector& alpha, sgpp::base::DataMatrix& data, sgpp::base::DataVector& classes,
    sgpp::base::DataVector& charaNumbers) {
  base::LinearBoundaryBasis<unsigned int, unsigned int> base;
  return test_datasetWithCharacteristicNumber(this->storage, base, alpha, data, classes,
                                              charaNumbers, 0);
}

void OperationTestLinearStretchedBoundary::calculateROCcurve(sgpp::base::DataVector& alpha,
                                                             sgpp::base::DataMatrix& data,
                                                             sgpp::base::DataVector& classes,
                                                             sgpp::base::DataVector& thresholds,
                                                             sgpp::base::DataMatrix& ROC_curve) {
  base::LinearBoundaryBasis<unsigned int, unsigned int> base;
  test_calculateROCcurve(this->storage, base, alpha, data, classes, thresholds, ROC_curve);
}
}  // namespace datadriven
}  // namespace sgpp
