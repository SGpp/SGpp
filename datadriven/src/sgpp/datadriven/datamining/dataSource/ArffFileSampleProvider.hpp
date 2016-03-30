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

#include <sgpp/datadriven/datamining/dataSource/SampleProvider.hpp>
#include <sgpp/datadriven/datamining/dataSource/FileSampleProvider.hpp>

#include <string>

namespace sgpp {
namespace datadriven {

class ArffFileSampleProvider : public FileSampleProvider {
 public:
  ArffFileSampleProvider();

  virtual ~ArffFileSampleProvider();

  std::unique_ptr<Dataset> getNextSamples(int howMany);

  std::unique_ptr<Dataset> getAllSamples();

  size_t getDim();

  size_t getDatasetSize();

  size_t getNumClasses();

  void readFile(std::string fileName);

  void readString(std::string& input);

 private:
  std::unique_ptr<Dataset> dataset;
};

} /* namespace datadriven */
} /* namespace sgpp */
