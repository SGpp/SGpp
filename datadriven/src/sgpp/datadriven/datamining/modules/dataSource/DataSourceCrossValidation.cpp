// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/DataSourceCrossValidation.hpp>

#include <vector>
#include <iostream>

using sgpp::base::DataVector;
using sgpp::base::algorithm_exception;

namespace sgpp {
namespace datadriven {

DataSourceCrossValidation::DataSourceCrossValidation(
    const DataSourceConfig& dataSourceConfig,
    const CrossvalidationConfiguration& crossValidationConfig,
    DataShufflingFunctorCrossValidation* shuffling, SampleProvider* sampleProvider)
    : DataSource{dataSourceConfig, sampleProvider},
      validationData{nullptr},
      crossValidationConfig{crossValidationConfig},
      shuffling(shuffling) {}

Dataset* DataSourceCrossValidation::getValidationData() { return validationData; }

void DataSourceCrossValidation::reset() {
  sampleProvider->reset();

  // Retrieve validation data again
  if (validationData != nullptr) delete validationData;
  size_t validationSize = shuffling->getCurrentFoldSize(sampleProvider->getNumSamples());
  validationData = sampleProvider->getNextSamples(validationSize);
}

void DataSourceCrossValidation::setFold(size_t foldIdx) { shuffling->setFold(foldIdx); }

const CrossvalidationConfiguration& DataSourceCrossValidation::getCrossValidationConfig() const {
  return crossValidationConfig;
}

} /* namespace datadriven */
} /* namespace sgpp */
