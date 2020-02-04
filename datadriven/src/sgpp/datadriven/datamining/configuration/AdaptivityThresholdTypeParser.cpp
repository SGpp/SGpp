// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/datadriven/datamining/configuration/AdaptivityThresholdTypeParser.hpp>

#include <algorithm>
#include <string>

namespace sgpp {
namespace datadriven {

using sgpp::base::data_exception;

AdaptivityThresholdType AdaptivityThresholdTypeParser::parse(const std::string& input) {
  auto inputLower = input;
  std::transform(inputLower.begin(), inputLower.end(), inputLower.begin(), ::tolower);
  if (inputLower == "relative") {
    return AdaptivityThresholdType::Relative;
  } else if (inputLower == "absolute") {
    return AdaptivityThresholdType::Absolute;
  } else {
    std::string errorMsg = "Failed to convert string \"" + input +
                           "\" to any known "
                           "AdaptivityThresholdType";
    throw data_exception(errorMsg.c_str());
  }
}

const std::string& AdaptivityThresholdTypeParser::toString(AdaptivityThresholdType type) {
  return refinementFunctorTypeMap.at(type);
}

const AdaptivityThresholdTypeParser::AdaptivityThresholdTypeMap_t
    AdaptivityThresholdTypeParser::refinementFunctorTypeMap = []() {
      return AdaptivityThresholdTypeMap_t{
          std::make_pair(AdaptivityThresholdType::Relative, "relative"),
          std::make_pair(AdaptivityThresholdType::Absolute, "absolute"),
      };
    }();
}  // namespace datadriven
}  // namespace sgpp
