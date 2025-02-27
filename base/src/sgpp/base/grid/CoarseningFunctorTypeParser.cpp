// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/CoarseningFunctorTypeParser.hpp>

#include <sgpp/base/exception/data_exception.hpp>
#include <algorithm>
#include <string>

namespace sgpp {
namespace base {

using sgpp::base::data_exception;

CoarseningFunctorType CoarseningFunctorTypeParser::parse(const std::string& input) {
  auto inputLower = input;
  std::transform(inputLower.begin(), inputLower.end(), inputLower.begin(), ::tolower);
  if (inputLower == "surplus") {
    return CoarseningFunctorType::Surplus;
  } else if (inputLower == "surplusvolume") {
    return CoarseningFunctorType::SurplusVolume;
  } else if (inputLower == "surplusabsolutevalue") {
    return CoarseningFunctorType::SurplusAbsoluteValue;
  } else if (inputLower == "classification") {
    return CoarseningFunctorType::Classification;
  } else {
    std::string errorMsg = "Failed to convert string \"" + input +
                           "\" to any known "
                           "CoarseningFunctorType";
    throw data_exception(errorMsg.c_str());
  }
}

const std::string& CoarseningFunctorTypeParser::toString(
    CoarseningFunctorType type) {
  return coarseningFunctorTypeMap.at(type);
}

const CoarseningFunctorTypeParser::CoarseningFunctorTypeMap_t
    CoarseningFunctorTypeParser::coarseningFunctorTypeMap = []() {
      return CoarseningFunctorTypeMap_t{
          std::make_pair(CoarseningFunctorType::Surplus, "Surplus"),
          std::make_pair(CoarseningFunctorType::SurplusVolume, "SurplusVolume"),
          std::make_pair(CoarseningFunctorType::SurplusAbsoluteValue, "SurplusAbsoluteValue"),
          std::make_pair(CoarseningFunctorType::Classification, "Classification")};
    }();

} /* namespace base */
} /* namespace sgpp */
