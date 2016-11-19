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
 * A FileSampleProvider generalizes the
 */
class FileSampleProvider : public SampleProvider {
 public:
  virtual ~FileSampleProvider() {}

  virtual void readFile(const std::string& fileName) = 0;
  virtual void readString(const std::string& input) = 0;
};

} /* namespace datadriven */
} /* namespace sgpp */
