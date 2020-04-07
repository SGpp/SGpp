// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/datadriven/configuration/CrossvalidationConfiguration.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/DataSource.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/shuffling/DataShufflingFunctorCrossValidation.hpp>

#include <vector>

namespace sgpp {
namespace datadriven {
/**
 * DataSourceCrossValidation is a high level interface to provide functionality
 * for processing
 * data using a cross validation enviroment. That is retrieving a certain fold
 * for validation
 * and the rest of the data for training.
 * Note that memory-wise this is very costly and not tractable for large data.
 */
class DataSourceCrossValidation : public DataSource {
 public:
  /**
   * Constructor
   * @param dataSourceConfig configuration of the data source
   * @param crossValidationConfig configuration of the cross validation
   * @param shuffling cross validation shuffling that is used by the sample provider instance
   * @param sampleProvider the sample provider to operate on.
   */
  DataSourceCrossValidation(const DataSourceConfig& dataSourceConfig,
                            const CrossvalidationConfiguration& crossValidationConfig,
                            DataShufflingFunctorCrossValidation* shuffling,
                            SampleProvider* sampleProvider);

  /**
   * Returns the data that is used for validation, i.e. the current fold.d If
   * all folds were already
   * iterated over, this method throws.
   * @return pointer to the validation dataset
   */
  Dataset* getValidationData() override;

  /**
   * Sets the next fold idx to be used for cross validation
   * @param foldIdx index of the fold
   */
  void setFold(size_t foldIdx);

  /**
   * Resets the state of the the sample provider to begin a new training epoch
   */
  void reset();

  /**
   * Gets the configuration for the cross validation.
   * @return configuration for the cross validation
   */
  const CrossvalidationConfiguration& getCrossValidationConfig() const;

  /**
   * Clean up memory
   */
  ~DataSourceCrossValidation() override {
    if (validationData != nullptr) {
      delete validationData;
    }
  }

 private:
  /**
   * Validation dataset
   */
  Dataset* validationData;
  /**
   * Configuration for the cross validation
   */
  CrossvalidationConfiguration crossValidationConfig;
  /**
   * Shuffling functor that is held by the sample provider.
   */
  DataShufflingFunctorCrossValidation* shuffling;
};

} /* namespace datadriven */
} /* namespace sgpp */
