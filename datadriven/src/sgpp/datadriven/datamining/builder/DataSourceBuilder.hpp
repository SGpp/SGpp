// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <memory.h>
#include <sgpp/datadriven/datamining/modules/dataSource/DataSource.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/DataSourceSplitting.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/DataSourceCrossValidation.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/DataSourceConfig.hpp>
#include <sgpp/datadriven/configuration/CrossvalidationConfiguration.hpp>
#include <string>
#include <vector>

namespace sgpp {
namespace datadriven {

/**
 * Generate an instance of #sgpp::datadriven::DataSource using the Builder Pattern.
 */
class DataSourceBuilder {
 public:
  /**
   * Default constructor
   */
  DataSourceBuilder() = default;

  /**
   * Optionally specify a valid path to a file that should be read if files are used.
   * @param filePath valid path to a file that should be read by the data source.
   * @return Reference to this object, used for chaining.
   */
  DataSourceBuilder& withPath(const std::string& filePath);

  /**
   * Optionally Specify if the file used is gz compressed. If data source does not use any files,
   * this is set to false by default.
   * @param isCompressed true if the file is compressed, false otherwise.
   * @return Reference to this object, used for chaining.
   */
  DataSourceBuilder& withCompression(bool isCompressed);

  /**
   * Optionally Specify the file type if files are used. If data source does not use any files,
   * this is set to none by default.
   * @param fileType value of
   * @return Reference to this object, used for chaining.
   */
  DataSourceBuilder& withFileType(DataSourceFileType fileType);

  /**
   * Optionally Specify the amount of batches if batch learning is used. If no batch learning is
   * used, all data is returned as a single batch (same as howMany=1).
   * @param howMany amount of batches used in batch learning scenario.
   * @return Reference to this object, used for chaining.
   */
  DataSourceBuilder& inBatches(size_t howMany);

  /**
   * Optionally Specify the batch size if batch learning is used. If no batch learning is used this
   * value defaults to 0 (all samples).
   * @param batchSize size of batches used in batch learning scenario.
   * @return Reference to this object, used for chaining.
   */
  DataSourceBuilder& withBatchSize(size_t batchSize);

  /**
   * Based on the currently specified configuration, build and configure an instance of a data
   * source object.
   * @return Fully configured instance of #sgpp::datadriven::DataSourceSplitting object.
   */
  DataSourceSplitting* splittingAssemble() const;

  /**
   * Factory method used to build an instance of a #sgpp::datadriven::DataSourceSplitting object
   * based on the passed configuration.
   * @param config configuration for the data source instance
   * @return Fully configured instance of #sgpp::datadriven::DataSourceSplitting object.
   */
  DataSourceSplitting* splittingFromConfig(const DataSourceConfig& config);

  /**
   * Based on the currently specified configuration, build and configure an instance of a data
   * source object that is able to perform cross validation.
   * @return Fully configured instance of #sgpp::datadriven::DataSourceCrossValidation object.
   */
  DataSourceCrossValidation* crossValidationAssemble() const;

  /**
   * Factory method used to build an instance of a #sgpp::datadriven::DataSourceCrossValidation
   * object based on the passed configuration.
   * @param config configuration for the data source instance
   * @param crossValidationConfig configuration for the cross validation
   * @return Fully configured instance of #sgpp::datadriven::DataSourceCrossValidation object.
   */
  DataSourceCrossValidation* crossValidationFromConfig(const DataSourceConfig& config,
      const CrossvalidationConfiguration &crossValidationConfig);

 private:
  /**
   * Extract file type and compression based on the file extensions.
   */
  void grabTypeInfoFromFilePath();

  /**
   * Current state of the object is stored inside this configuration object.
   */
  DataSourceConfig config;

  /**
   * Current state of the object is stored inside this configuration object (for cross validation).
   */
  CrossvalidationConfiguration crossValidationConfig;
};

} /* namespace datadriven */
} /* namespace sgpp */
