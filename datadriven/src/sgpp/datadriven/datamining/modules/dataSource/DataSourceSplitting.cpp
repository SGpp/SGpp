/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * DataSourceSplitting.cpp
 *
 *  Created on: Jul 23, 2018
 *      Author: dominik
 */

#include <sgpp/datadriven/datamining/modules/dataSource/DataSourceSplitting.hpp>

namespace sgpp {
namespace datadriven {

DataSourceSplitting::DataSourceSplitting(DataSourceConfig config, SampleProvider* sampleProvider) :
  DataSource{config, sampleProvider} {
  size_t validationSize = static_cast<size_t>(config.validationPortion *
      static_cast<double>(sampleProvider->getNumSamples()));
  validationData = sampleProvider->getNextSamples(validationSize);
}

Dataset *DataSourceSplitting::getValidationData() { return validationData; }

void DataSourceSplitting::reset() {
  sampleProvider->reset();
}

} /* namespace datadriven */
} /* namespace sgpp */




