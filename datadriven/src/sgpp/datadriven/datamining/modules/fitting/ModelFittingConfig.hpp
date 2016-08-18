/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * ModelFittingConfig.hpp
 *
 * Created on: Aug 18, 2016
 *     Author: Michael Lettrich
 */

#pragma once

#include <sgpp/base/exception/data_exception.hpp>

using sgpp::base::data_exception;

namespace sgpp {
namespace datadriven {

enum class FitterType { RegressionLeastSquares };

class FitterTypeParser {
 public:
  static FitterType parse(const std::string& input) {
    auto inputLower = input;
    std::transform(inputLower.begin(), inputLower.end(), inputLower.begin(), ::tolower);

    if (inputLower == "regressionleastsquares") {
      return FitterType::RegressionLeastSquares;
    } else {
      std::string errorMsg = "Failed to convert string \"" + input + "\" to any known FitterType";
      throw data_exception(errorMsg.c_str());
    }
  };
};

} /* namespace datadriven */
} /* namespace sgpp */
