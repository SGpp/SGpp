// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/modules/dataSource/DataSourceSplitting.hpp>

namespace sgpp {
namespace datadriven {

DataSourceSplitting::DataSourceSplitting(const DataSourceConfig& config,
                                         SampleProvider* sampleProvider)
    : DataSource{config, sampleProvider}, validationData{nullptr} {}

Dataset* DataSourceSplitting::getValidationData() { return validationData; }

void DataSourceSplitting::reset() {
  sampleProvider->reset();
  // Retrieve new validation data
  delete validationData;
  size_t validationSize = static_cast<size_t>(config.validationPortion_ *
                                              static_cast<double>(sampleProvider->getNumSamples()));
  validationData = sampleProvider->getNextSamples(validationSize);
}

} /* namespace datadriven */
} /* namespace sgpp */
