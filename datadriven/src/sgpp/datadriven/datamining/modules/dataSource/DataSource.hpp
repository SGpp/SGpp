// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/datadriven/datamining/modules/dataSource/DataSourceConfig.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/DataSourceIterator.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/DataTransformation.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/SampleProvider.hpp>
#include <sgpp/datadriven/tools/Dataset.hpp>

#include <string>

namespace sgpp {
namespace datadriven {

// forward declaration
class DataSourceIterator;

/**
 * DataSource is a high level, easy to use interface for accessing data provided by a all kinds of
 * #sgpp::datadriven::SampleProvider. Should be used by end users.
 */
class DataSource {
 public:
  /**
   * Constructor
   * @param config configuration object used for the data source
   * @param sampleProvider the sample provider to operate on.
   */
  DataSource(DataSourceConfig config, SampleProvider* sampleProvider);

  virtual ~DataSource() = default;

  /**
   * Read only access to the configuration used by DataSource and underlying SampleProvider.
   * @return Current configuration object.
   */
  const DataSourceConfig& getConfig() const;

  /**
   * Request data from the underlying SampleProvider as specified in the provided configuration
   * object upon construction.
   * @return #sgpp::datadriven::Dataset containing requested amount of samples (if available).
   */
  virtual Dataset* getNextSamples();

  /**
   * Return an iterator object pointing to the first batch of this DataSource. Can be used to obtain
   * new batches in batch learning scenarios as often as specified inside the configuration. Allows
   * convenient range based for loops for batch learning.
   * @return iterator object pointing to the first batch.
   */
  DataSourceIterator begin();

  /**
   * Return an iterator object pointing to the last possible batch of this DataSource. Required for
   * range based for loops.
   * @return iterator object pointing to the last possible batch.
   */
  DataSourceIterator end();

  /**
   * Return how many batches have already been requested from this DataSource. Required for
   * range based for loops using the DataSourceIterator.
   * @return the amount of batches that have already been requested.
   */
  size_t getCurrentIteration() const;

  /**
   * Returns the data that is used for validation
   * @return pointer to the validation dataset
   */
  virtual Dataset *getValidationData() = 0;

 protected:
  /**
   * Configuration file that determines all relevant properties of the object.
   */
  DataSourceConfig config;

  /**
   * counter variable if data is requested in batches.
   */
  size_t currentIteration;

  /**
   * pointer to sample provider that actually handles data aquisition.
   */
  std::unique_ptr<SampleProvider> sampleProvider;

  /**
   * pointer to DataTransformation to perform transformations on init.
   */
  DataTransformation* dataTransformation;
};

} /* namespace datadriven */
} /* namespace sgpp */
