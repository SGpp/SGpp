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

#include <sgpp/base/exception/data_exception.hpp>

#include <algorithm>
#include <string>

using sgpp::base::data_exception;

namespace sgpp {
namespace datadriven {

enum class DataSourceFileType { NONE, ARFF };

class DataSourceFileTypeParser {
 public:
  DataSourceFileType operator()(const std::string& input) const {
    auto inputLower = input;
    std::transform(inputLower.begin(), inputLower.end(), inputLower.begin(), ::tolower);

    if (inputLower == arff) {
      return DataSourceFileType::ARFF;
    } else if (inputLower == none) {
      return DataSourceFileType::NONE;
    } else {
      std::string errorMsg =
          "Failed to convert string \"" + input + "\" to any known DataSourceFileType";
      throw data_exception(errorMsg.c_str());
    }
  };

 private:
  const std::string arff = "arff";
  const std::string none = "none";
};

struct DataSourceConfig {
  std::string filePath = "";
  DataSourceFileType fileType = DataSourceFileType::NONE;
  bool isCompressed = false;
  size_t numBatches = 1;
  size_t batchSize = 0;
};

} /* namespace datadriven */
} /* namespace sgpp */
