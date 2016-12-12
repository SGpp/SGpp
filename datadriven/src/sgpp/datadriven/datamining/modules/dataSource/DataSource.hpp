/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * SampleProviderModule.hpp
 *
 *  Created on: 17.05.2016
 *      Author: Michael Lettrich
 */

#pragma once

#include <sgpp/datadriven/datamining/modules/dataSource/DataSourceConfig.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/DataSourceIterator.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/SampleProvider.hpp>
#include <sgpp/datadriven/tools/Dataset.hpp>

#include <string>

namespace sgpp {
namespace datadriven {

// forward declaration
class DataSourceIterator;

/**
 * Configurable high level object that controlls
 */
class DataSource {
 public:
  /**
   *
   * @param config
   * @param sampleProvider
   */
  DataSource(DataSourceConfig config, SampleProvider* sampleProvider);

  /**
   *
   */
  const DataSourceConfig& getConfig() const;

  /**
   *
   */
  Dataset* getNextSamples();
  /**
   *
   */
  DataSourceIterator begin();
  /**
   *
   */
  DataSourceIterator end();

  /**
   *
   */
  size_t getCurrentIteration() const;
  size_t getCurrentIteration();

 private:
  /**
   * Configuration file that determines all relevant properties of the object.
   */
  DataSourceConfig config;

  /**
   * counter variable if data is requested in batches.
   */
  size_t currentIteration;

  /**
   * pointer to sample provider that actually returns the data.
   */
  std::unique_ptr<SampleProvider> sampleProvider;
};

} /* namespace datadriven */
} /* namespace sgpp */
