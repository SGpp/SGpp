/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * SampleProviderIterator.hpp
 *
 *  Created on: 17.05.2016
 *      Author: Michael Lettrich
 */

#pragma once

#include <sgpp/datadriven/datamining/modules/dataSource/DataSource.hpp>
#include <sgpp/datadriven/tools/Dataset.hpp>

namespace sgpp {
namespace datadriven {
// forward declaration to break dependency cycle
class DataSource;

/**
 * Iterator object for walking convenient walking over the batches of a
 * #sgpp::datadriven::Datasource
 */
class DataSourceIterator {
 public:
  /**
   *
   */
  DataSourceIterator(DataSource& sampleProvider, size_t counter);

  /**
   *
   */
  bool operator!=(const DataSourceIterator& other);

  /**
   *
   */
  const DataSourceIterator& operator++();

  /**
   *
   */
  Dataset* operator*() const;

 private:
  /**
   * The sample provider to iterate over
   */
  DataSource& sampleProvider;

  /**
   * Private counter to keep track of the current iteration.
   */
  size_t counter;
};

} /* namespace datadriven */
} /* namespace sgpp */
