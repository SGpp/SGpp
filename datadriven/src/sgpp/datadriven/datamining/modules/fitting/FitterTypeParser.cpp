// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/modules/fitting/FitterTypeParser.hpp>

#include <sgpp/base/exception/data_exception.hpp>

#include <algorithm>
#include <string>

namespace sgpp {
namespace datadriven {

using sgpp::base::data_exception;

FitterType FitterTypeParser::parse(const std::string &input) {
  auto inputLower = input;
  std::transform(inputLower.begin(), inputLower.end(), inputLower.begin(), ::tolower);

  if (inputLower == "regressionleastsquares") {
    return FitterType::RegressionLeastSquares;
  } else if (inputLower == "densityestimation") {
    return FitterType::DensityEstimation;
  } else if (inputLower == "densityratioestimation") {
    return FitterType::DensityRatioEstimation;
  } else if (inputLower == "densitydifferenceestimation") {
    return FitterType::DensityDifferenceEstimation;
  } else if (inputLower == "classification") {
    return FitterType::Classification;
  } else {
    std::string errorMsg = "Failed to convert string \"" + input + "\" to any known FitterType";
    throw data_exception(errorMsg.c_str());
  }
}

const std::string &FitterTypeParser::toString(FitterType type) { return fitterTypeMap.at(type); }

const FitterTypeParser::FitterTypeMap_t FitterTypeParser::fitterTypeMap = []() {
  return FitterTypeParser::FitterTypeMap_t{
      std::make_pair(FitterType::RegressionLeastSquares, "ModelFittingLeastSquares"),
      std::make_pair(FitterType::DensityEstimation, "ModelFittingDensityEstimation"),
      std::make_pair(FitterType::DensityRatioEstimation, "ModelFittingDensityRatioEstimation"),
      std::make_pair(FitterType::DensityDifferenceEstimation,
                     "ModelFittingDensityDifferenceEstimation"),
      std::make_pair(FitterType::Classification, "ModelFittingClassification")};
}();

} /* namespace datadriven */
} /* namespace sgpp */
