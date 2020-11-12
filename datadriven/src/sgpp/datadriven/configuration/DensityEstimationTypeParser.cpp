// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/configuration/DensityEstimationTypeParser.hpp>

#include <sgpp/base/exception/data_exception.hpp>
#include <algorithm>
#include <string>

namespace sgpp {
namespace datadriven {

using sgpp::base::data_exception;

DensityEstimationType DensityEstimationTypeParser::parse(
    const std::string &input) {
  auto inputLower = input;
  std::transform(inputLower.begin(), inputLower.end(), inputLower.begin(),
                 ::tolower);

  if (inputLower.compare("cg") == 0) {
    return sgpp::datadriven::DensityEstimationType::CG;
  } else if (inputLower.compare("decomposition") == 0) {
    return sgpp::datadriven::DensityEstimationType::Decomposition;
  } else {
    std::string errorMsg = "Failed to convert string \"" + input +
                           "\" to any known "
                           "DensityEstimationType";
    throw data_exception(errorMsg.c_str());
  }
}

const std::string &DensityEstimationTypeParser::toString(
    DensityEstimationType type) {
  return densityEstimationTypeMap.at(type);
}

const DensityEstimationTypeParser::DensityEstimationTypeMap_t
    DensityEstimationTypeParser::densityEstimationTypeMap = []() {
      return DensityEstimationTypeMap_t{
          std::make_pair(DensityEstimationType::CG, "CG"),
          std::make_pair(DensityEstimationType::Decomposition,
                         "Decomposition")};
    }();
} /* namespace datadriven */
} /* namespace sgpp */
