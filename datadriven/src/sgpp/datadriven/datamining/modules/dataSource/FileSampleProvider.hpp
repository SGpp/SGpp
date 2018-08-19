/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * FileSampleProvider.hpp
 *
 *  Created on: 10.03.2016
 *      Author: Michael Lettrich
 */

#pragma once

#include <sgpp/datadriven/datamining/modules/dataSource/SampleProvider.hpp>

#include <string>

namespace sgpp {
namespace datadriven {

/**
 * #sgpp::datadriven::FileSampleProvider is an specialization of #sgpp::datadriven::SampleProvider
 * and provides an interface for all sample providers that get their samples from files.
 */
class FileSampleProvider : public SampleProvider {
 public:
  /**
   * Read the contents of the file at the given path. Has to throw an exception if file can not be
   * opened or parsed. Results of parsing can be optained via #sgpp::datadriven::SampleProvider
   * member functions.
   * @param filePath valid path to an existing file.
   * @param hasTargets whether the file has targets (i.e. supervised learning)
   */
  virtual void readFile(const std::string &filePath, bool hasTargets) = 0;

  /**
   * Read the contents of a string, for example a deflated archive. Has to throw an exception if
   * string can not be parsed. Results of parsing can be optained via
   * #sgpp::datadriven::SampleProvider member functions.
   * @param input the raw string input to parse
   * @param hasTargets whether the file has targest (i.e. supervised learning)
   */
  virtual void readString(const std::string &input, bool hasTargets) = 0;
};
} /* namespace datadriven */
} /* namespace sgpp */
