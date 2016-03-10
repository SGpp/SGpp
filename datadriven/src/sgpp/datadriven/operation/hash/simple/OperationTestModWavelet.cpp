// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/common/basis/PolyModifiedBasis.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationTestModWavelet.hpp>

#include <sgpp/datadriven/algorithm/test_dataset.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {

double OperationTestModWavelet::test(sgpp::base::DataVector& alpha, sgpp::base::DataMatrix& data,
                                      sgpp::base::DataVector& classes) {
  base::WaveletModifiedBasis<unsigned int, unsigned int> base;
  return test_dataset(this->storage, base, alpha, data, classes);
}

double OperationTestModWavelet::testMSE(sgpp::base::DataVector& alpha,
                                         sgpp::base::DataMatrix& data,
                                         sgpp::base::DataVector& refValues) {
  base::WaveletModifiedBasis<unsigned int, unsigned int> base;
  return test_dataset_mse(this->storage, base, alpha, data, refValues);
}

double OperationTestModWavelet::testWithCharacteristicNumber(
    sgpp::base::DataVector& alpha, sgpp::base::DataMatrix& data, sgpp::base::DataVector& classes,
    sgpp::base::DataVector& charaNumbers) {
  base::WaveletModifiedBasis<unsigned int, unsigned int> base;
  return test_datasetWithCharacteristicNumber(this->storage, base, alpha, data, classes,
                                              charaNumbers, 0.0);
}

void OperationTestModWavelet::calculateROCcurve(sgpp::base::DataVector& alpha,
                                                sgpp::base::DataMatrix& data,
                                                sgpp::base::DataVector& classes,
                                                sgpp::base::DataVector& thresholds,
                                                sgpp::base::DataMatrix& ROC_curve) {
  base::WaveletModifiedBasis<unsigned int, unsigned int> base;
  test_calculateROCcurve(this->storage, base, alpha, data, classes, thresholds, ROC_curve);
}

}  // namespace datadriven
}  // namespace sgpp
