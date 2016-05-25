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

#include <sgpp/datadriven/datamining/modules/dataSource/FileSampleProvider.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/SampleProviderIterator.hpp>
#include <sgpp/datadriven/tools/Dataset.hpp>
#include <sgpp/globaldef.hpp>

#include <limits>
#include <memory>
#include "DataSource.hpp"

namespace sgpp {
namespace datadriven {

DataSource::DataSource(std::shared_ptr<DataSourceState> state,
                       std::unique_ptr<SampleProvider> sampleProvider)
    : state(state), sampleProvider(std::move(sampleProvider)) {
  // if a file name was specified, we are reading from a file, so we need to open it.
  if (!state->getFilePath().empty()) {
    FileSampleProvider* fileProvider = dynamic_cast<FileSampleProvider*>(sampleProvider.get());
    fileProvider->readFile(state->getFilePath());
  }
}

DataSource::~DataSource() {}

DataSourceIterator DataSource::begin() { return DataSourceIterator(*this, 0); }

DataSourceIterator DataSource::end() { return DataSourceIterator(*this, state->getNumBatches()); }

std::unique_ptr<Dataset> DataSource::getNextSamples() {
  // only one iteration: we want all samples
  if (state->getNumBatches() == 1 && state->getBatchSize() != 0) {
    state->incrementCurrentIteration();
    return sampleProvider->getAllSamples();
  } else {
    state->incrementCurrentIteration();
    return sampleProvider->getNextSamples(state->getBatchSize());
  }
}

} /* namespace datadriven */
} /* namespace sgpp */
