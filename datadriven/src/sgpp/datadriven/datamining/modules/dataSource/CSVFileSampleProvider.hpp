/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * CSVFileSampleProvider.hpp
 *
 *  Created on: Dec 18, 2017
 *      Author: Eric Koepke, perun, Michael Lettrich
 */

#pragma once

#include <sgpp/datadriven/datamining/modules/dataSource/FileSampleProvider.hpp>

#include <string>

// TODO(lettrich): allow different splitting techniques e.g. proportional splitting for
// classification
// TODO(lettrich): save memory by only storing the dataset if data is accessed in batches.
namespace sgpp {
namespace datadriven {

/**
 * CSVFileSampleProvider allows reading data in CSV format into a #sgpp::datadriven::Dataset
 * object. Data can currently only be a file containing CSV data with the first line containing
 * column titles (is skipped).
 *
 */
class CSVFileSampleProvider : public FileSampleProvider {
 public:
  /**
   * Default constructor
   * @param validationPortion portion of data that is used for validation
   */
  explicit CSVFileSampleProvider(double validationPortion);

  /**
   * Clone Pattern to allow copying of derived classes.
   * @return a Pointer to a new instance of #sgpp::datadriven::CSVFileSampleProvider with copied
   * state. Caller owns the new object.
   */
  SampleProvider* clone() const override;

  Dataset* getNextSamples(size_t howMany) override;

  Dataset* getAllSamples() override;

  size_t getDim() const override;

  size_t getNumSamples() const override;

  /**
   * Open an existing CSV file, parse it and store its contents inside this class. Throws if file
   * can not be opened or parsed.
   * @param filePath Path to an existing file.
   * @param hasTargets whether the file has targest (i.e. supervised learning)
   */
  void readFile(const std::string& filePath, bool hasTargets) override;

  /**
   * Currently not implemented.
   * @param input string containing information in CSV file format
   * @param hasTargets whether the file has targest (i.e. supervised learning)
   */
  void readString(const std::string& input, bool hasTargets) override;

  /**
   * Returns the data that is used for validation
   * @return pointer to the validation dataset
   */
  Dataset *getValidationData();

 private:
  /**
   * Portion of the dataset that is used for validation
   */
  double validationPortion;

  /**
   * #sgpp::datadriven::Dataset containing the samples read from file or string.
   */
  Dataset dataset;

  /**
   * Indicates the index in dataset where #getNextSamples will start grabbing new samples in its
   * next call. After each call of #getNextSamples, the counter is set to the amount of min(counter
   * + requestedSamplesSize, dataset.getSize()).
   */
  size_t counter;

  /**
   * Helper member function for #getNextSamples. Linearly walks through dataset, beginning at
   * counter and returns a pointer to a new instance of #sgpp::datadriven::Dataset containing the
   * desired amount of samples (if available - else all remaining samples) and updates counter.
   */
  Dataset* splitDataset(size_t howMany);
};
} /* namespace datadriven */
} /* namespace sgpp */
