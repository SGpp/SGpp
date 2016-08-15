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

#include <string>

namespace sgpp {
namespace datadriven {

enum class DataSourceFileType { NONE, ARFF };

struct DataSourceConfig {
  std::string filePath = "";
  DataSourceFileType fileType = DataSourceFileType::NONE;
  bool isCompressed = false;
  size_t numBatches = 1;
  size_t batchSize = 0;
};
} /* namespace datadriven */
} /* namespace sgpp */
