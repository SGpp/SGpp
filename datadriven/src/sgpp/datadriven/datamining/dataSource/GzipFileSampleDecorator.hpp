/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * GzipFileSampleDecorator.hpp
 *
 *  Created on: 01.04.2016
 *      Author: Michael Lettrich
 */

#pragma once

#include <sgpp/datadriven/datamining/dataSource/FileSampleDecorator.hpp>
#include <memory>

namespace sgpp {
namespace datadriven {

class GzipFileSampleDecorator : public FileSampleDecorator {
 public:
  explicit GzipFileSampleDecorator(std::unique_ptr<FileSampleProvider> fileSampleProvider);

  virtual ~GzipFileSampleDecorator();

  std::unique_ptr<Dataset> getNextSamples(size_t howMany);

  std::unique_ptr<Dataset> getAllSamples();

  size_t getDim();

  size_t getDatasetSize();

  // size_t getNumClasses();

  void readFile(const std::string& fileName);

  void readString(const std::string& input);
};

} /* namespace datadriven */
} /* namespace sgpp */
