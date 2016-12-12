/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * ARFFWrapper.hpp
 *
 *  Created on: Feb 8, 2016
 *      Author: perun, Michael Lettrich
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
 * ArffFileSampleProvider allows reading data in ARFF format into a #sgpp::datadriven::Dataset
 * object. Data can currently be either be a string formatted in ARFF or a file containing ARFF
 * data.
 *
 */
class ArffFileSampleProvider : public FileSampleProvider {
 public:
  /**
   * Default constructor
   */
  ArffFileSampleProvider();

  /**
   * Clone Pattern to allow copying of derived classes.
   * @return a Pointer to a new instance of #sgpp::datadriven::ArffFileSampleProvider with copied
   * state. Caller owns the new object.
   */
  SampleProvider* clone() const override;

  Dataset* getNextSamples(size_t howMany) override;

  Dataset* getAllSamples() override;

  size_t getDim() const override;

  size_t getNumSamples() const override;

  /**
   * Open an existing ARFF file, parse it and store its contents inside this class. Throws if file
   * can not be opened or parsed.
   * @param filePath Path to an existing file.
   */
  void readFile(const std::string& filePath) override;

  /**
   * Parse contents of a string containing information in ARFF format, parse it and store its
   * contents inside this class. Throws if string can not be parsed.
   * @param input string containing information in ARFF file format
   */
  void readString(const std::string& input) override;

 private:
  /**
   * #sgpp::datadriven::Dataset containing the samples read from file or string.
   */
  Dataset dataset;

  /**
   * Indicates the index in #dataset where #getNextSamples will start grabing new samples in its
   * next call. After each call of #getNextSamples, the counter is set to the amount of min(#counter
   * + requestedSamplesSize, #dataset.getSize()).
   */
  size_t counter;

  /**
   * Helper member function for #getNextSamples. Linearly walks through #dataset, beginning at
   * #counter and returns a pointer to a new instance of #sgpp::datadriven::Dataset containing the
   * desired amount of samples (if available - else all remaining samples) and updates #counter.
   */
  Dataset* splitDataset(size_t howMany);
};

} /* namespace datadriven */
} /* namespace sgpp */
