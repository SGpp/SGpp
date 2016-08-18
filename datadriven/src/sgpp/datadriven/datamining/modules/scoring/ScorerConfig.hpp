/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * ScorerConfig.hpp
 *
 * Created on: Aug 18, 2016
 *     Author: Michael Lettrich
 */

#pragma once

#include <sgpp/base/exception/data_exception.hpp>

using sgpp::base::data_exception;

namespace sgpp {
namespace datadriven {

enum class ScorerShufflingType { random, sequential };
enum class ScorerMetric { MSE };

struct CrossValidationConfiguration {
  size_t folds = 5;
  int64_t randomSeed = -1;
  ScorerMetric metric;
  ScorerShufflingType shuffling = ScorerShufflingType::random;
};

struct TestingConfiguration {
  double testingPortion = 0.3;
  int64_t randomSeed = -1;
  ScorerMetric metric;
  ScorerShufflingType shuffling = ScorerShufflingType::random;
};

class ScorerShufflingTypeParser {
 public:
  static ScorerShufflingType parse(const std::string& input) {
    auto inputLower = input;
    std::transform(inputLower.begin(), inputLower.end(), inputLower.begin(), ::tolower);

    if (inputLower == "random") {
      return ScorerShufflingType::random;
    } else if (inputLower == "sequential") {
      return ScorerShufflingType::sequential;
    } else {
      std::string errorMsg =
          "Failed to convert string \"" + input + "\" to any known ScorerShufflingType";
      throw data_exception(errorMsg.c_str());
    }
  };
};

class ScorerMetricParser {
 public:
  static ScorerMetric parse(const std::string& input) {
    auto inputLower = input;
    std::transform(inputLower.begin(), inputLower.end(), inputLower.begin(), ::tolower);

    if (inputLower == "mse") {
      return ScorerMetric::MSE;
    } else {
      std::string errorMsg = "Failed to convert string \"" + input + "\" to any known ScorerMetric";
      throw data_exception(errorMsg.c_str());
    }
  };
};

} /* namespace datadriven */
} /* namespace sgpp */
