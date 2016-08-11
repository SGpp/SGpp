/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * SampleProviderModule.cpp
 *
 *  Created on: 17.05.2016
 *      Author: Michael Lettrich
 */

#include "DataSource.hpp"

#include <sgpp/datadriven/datamining/modules/dataSource/DataSourceIterator.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/FileSampleProvider.hpp>
#include <sgpp/datadriven/tools/Dataset.hpp>
#include <sgpp/globaldef.hpp>

#include <iostream>
#include <limits>
#include <memory>

namespace sgpp {
namespace datadriven {

DataSource::DataSource(DataSourceConfig conf, SampleProvider* sp)
    : config(conf), sampleProvider(std::shared_ptr<SampleProvider>(sp)) {
  // if a file name was specified, we are reading from a file, so we need to open it.
  if (!this->config.getFilePath().empty()) {
    auto fileProvider = std::static_pointer_cast<FileSampleProvider>(sampleProvider);
    fileProvider->readFile(this->config.getFilePath());
  }
}

DataSource::~DataSource() {}

DataSourceIterator DataSource::begin() { return DataSourceIterator(*this, 0); }

DataSourceIterator DataSource::end() { return DataSourceIterator(*this, config.getNumBatches()); }

Dataset* DataSource::getNextSamples() {
  // only one iteration: we want all samples
  if (config.getNumBatches() == 1 && config.getBatchSize() == 0) {
    config.incrementCurrentIteration();
    return sampleProvider->getAllSamples();
  } else {
    config.incrementCurrentIteration();
    return sampleProvider->getNextSamples(config.getBatchSize());
  }
}

DataSourceConfig& sgpp::datadriven::DataSource::getConfig() { return config; }

} /* namespace datadriven */
} /* namespace sgpp */
