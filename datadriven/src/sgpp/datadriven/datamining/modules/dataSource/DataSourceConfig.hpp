/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * DataSourceState.hpp
 *
 *  Created on: 25.05.2016
 *      Author: Michael Lettrich
 */

#pragma once

#include <sgpp/datadriven/datamining/modules/dataSource/DataTransformationConfig.hpp>

#include <string>

namespace sgpp {
namespace datadriven {

/**
 * Supported file types for sgpp::datadriven::FileSampleProvider
 */
enum class DataSourceFileType { NONE, ARFF, CSV };

/**
 * Enumeration of all supported shuffling types used to permute samples in a dataset. An entry
 * exists for each object that derives from #sgpp::datadriven::DataShufflingFunctor. Used for
 * configuration and factory methods.
 */
enum class DataSourceShufflingType { random, sequential };

/**
 * Configuration structure used for all kinds of SampleProviders including default values.
 */
struct DataSourceConfig {
  /**
   * Valid path to a file on disk. Empty for generated artificial datasets
   */
  std::string filePath = "";
  /**
   * Which type of input file are we dealing with? NONE for auto detection or generated artificial
   * datasets.
   */
  DataSourceFileType fileType = DataSourceFileType::NONE;
  /**
   * The dataset is gzip compressed
   */
  bool isCompressed = false;
  /**
   * How many batches should the dataset be split into for batch learning - if 1, take the
   * entire dataset
   */
  size_t numBatches = 1;
  /*
   * size of a batch - if 0, take all available samples.
   */
  size_t batchSize = 0;
  /*
   * The portion of the dataset that is used for validation
   */
  double validationPortion = 0.3;
  /**
   * whether the file has targets (i.e. supervised learning)
   */
  bool hasTargets = true;
  /*
   * Configuration for possible data transformation on dataset
   */
  datadriven::DataTransformationConfig dataTransformationConfig;
  /**
   * The type of shuffling to be applied to the data
   */
  DataSourceShufflingType shuffling = DataSourceShufflingType::sequential;
  /**
   * Seed for the shuffling prng
   */
  int64_t randomSeed = -1;
  /**
   * The number of epochs to train on
   */
  size_t epochs = 1;
};
} /* namespace datadriven */
} /* namespace sgpp */
