/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * RefinementFunctorTypeParser.cpp
 *
 *  Created on: Jun 29, 2018
 *
 */

#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/datadriven/datamining/configuration/RefinementFunctorTypeParser.hpp>

#include <algorithm>
#include <string>

namespace sgpp {
namespace datadriven {

using sgpp::base::data_exception;

RefinementFunctorType RefinementFunctorTypeParser::parse(const std::string& input) {
  auto inputLower = input;
  std::transform(inputLower.begin(), inputLower.end(), inputLower.begin(), ::tolower);
  if (inputLower == "surplus") {
    return RefinementFunctorType::Surplus;
  } else if (inputLower == "surplusvolume") {
    return RefinementFunctorType::SurplusVolume;
  } else if (inputLower == "zerocrossing") {
    return RefinementFunctorType::ZeroCrossing;
  } else if (inputLower == "databased") {
    return RefinementFunctorType::DataBased;
  } else {
    std::string errorMsg = "Failed to convert string \"" + input + "\" to any known "
        "RefinementFunctorType";
    throw data_exception(errorMsg.c_str());
  }
}

const std::string& RefinementFunctorTypeParser::toString(RefinementFunctorType type)
{ return refinementFunctorTypeMap.at(type); }

const RefinementFunctorTypeParser::RefinementFunctorTypeMap_t
RefinementFunctorTypeParser::refinementFunctorTypeMap = []() {
  return RefinementFunctorTypeMap_t{
      std::make_pair(RefinementFunctorType::Surplus, "Surplus"),
      std::make_pair(RefinementFunctorType::SurplusVolume, "SurplusVolume"),
      std::make_pair(RefinementFunctorType::ZeroCrossing, "ZeroCrossing"),
      std::make_pair(RefinementFunctorType::DataBased, "DataBased")
  };
}();

} /* namespace datadriven */
} /* namespace sgpp */




