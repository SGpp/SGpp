// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/ScorerMetricTypeParser.hpp>

#include <algorithm>
#include <map>
#include <string>

namespace sgpp {
namespace datadriven {

using sgpp::base::data_exception;

const ScorerMetricTypeParser::MetricTypeMap_t ScorerMetricTypeParser::metricTypeMap = []() {
  return MetricTypeMap_t{std::make_pair(ScorerMetricType::mse, "MSE"),
  std::make_pair(ScorerMetricType::nll, "NLL"),
  std::make_pair(ScorerMetricType::accuracy, "Accuracy")};
}();

const std::string& ScorerMetricTypeParser::toString(ScorerMetricType type) {
  return metricTypeMap.at(type);
}
ScorerMetricType ScorerMetricTypeParser::parse(const std::string& input) {
  auto inputLower = input;
  std::transform(inputLower.begin(), inputLower.end(), inputLower.begin(), ::tolower);

  if (inputLower == "mse") {
    return ScorerMetricType::mse;
  } else if (inputLower == "nll") {
    return ScorerMetricType::nll;
  } else if (inputLower == "accuracy") {
    return ScorerMetricType::accuracy;
  } else {
    const auto errorMsg = "Failed to convert string \"" + input + "\" to any known ScorerMetric";
    throw data_exception(errorMsg.c_str());
  }
}

} /* namespace datadriven */
} /* namespace sgpp */
