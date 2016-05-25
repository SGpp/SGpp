/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * SampleProviderIterator.cpp
 *
 *  Created on: 17.05.2016
 *      Author: Michael Lettrich
 */

#include "DataSourceIterator.hpp"

namespace sgpp {
namespace datadriven {
DataSourceIterator::DataSourceIterator(DataSource& sampleProvider, size_t counter)
    : sampleProvider(sampleProvider), counter(counter) {}

bool DataSourceIterator::operator!=(const DataSourceIterator& other) {
  return (counter != other.counter);
}

const DataSourceIterator& DataSourceIterator::operator++() {
  counter++;
  return *this;
}

std::unique_ptr<Dataset> DataSourceIterator::operator*() const {
  return sampleProvider.getNextSamples();
}

} /* namespace datadriven */
} /* namespace sgpp */
