/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * CoarseningFunctorTypeParser.cpp
 *
 *  Created on: May 28, 2019
 *
 */
#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/datadriven/datamining/configuration/CoarseningFunctorTypeParser.hpp>

#include <algorithm>
#include <string>

namespace sgpp {
namespace datadriven {

using sgpp::base::data_exception;

CoarseningFunctorType CoarseningFunctorTypeParser::parse(const std::string& input) {
  auto inputLower = input;
  std::transform(inputLower.begin(), inputLower.end(), inputLower.begin(), ::tolower);
  if (inputLower == "gridpointbased") {
    return CoarseningFunctorType::GridPointBased;
  } else {
    std::string errorMsg = "Failed to convert string \"" + input + "\" to any known "
        "CoarseningFunctorType";
    throw data_exception(errorMsg.c_str());
  }
}

const std::string& CoarseningFunctorTypeParser::toString(CoarseningFunctorType type)
{ return coarseningFunctorTypeMap.at(type); }

const CoarseningFunctorTypeParser::CoarseningFunctorTypeMap_t
CoarseningFunctorTypeParser::coarseningFunctorTypeMap = []() {
  return CoarseningFunctorTypeMap_t{
      std::make_pair(CoarseningFunctorType::GridPointBased, "GridPointBased")
  };
}();

} /* namespace datadriven */
} /* namespace sgpp */




