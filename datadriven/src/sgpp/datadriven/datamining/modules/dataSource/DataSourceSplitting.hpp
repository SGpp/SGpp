// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/datadriven/datamining/modules/dataSource/DataSource.hpp>

namespace sgpp {
namespace datadriven {
/**
 * DataSourceSlitting is a high level interface to provide functionality for
 * processing data
 * epoch-wise with a validation set that is retrieved at initialization time
 * using the first
 * samples the sample provider provides.
 */
class DataSourceSplitting : public DataSource {
 public:
  /**
   * Constructor
   * @param config configuration object used for the data source
   * @param sampleProvider the sample provider to operate on.
   */
  DataSourceSplitting(const DataSourceConfig& config, SampleProvider* sampleProvider);

  /**
   * Returns the data that is used for validation
   * @return pointer to the validation dataset
   */
  Dataset* getValidationData() override;

  /**
   * Resets the state of the the sample provider to begin a new training epoch
   */
  void reset();

  /**
   * Clean up memory
   */
  ~DataSourceSplitting() override {
    if (validationData != nullptr) {
      delete validationData;
    }
  }

 private:
  /**
   * The validation data that is retrieved from the sample provider at
   * initialization time
   */
  Dataset* validationData;
};

} /* namespace datadriven */
} /* namespace sgpp */
