/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * GeometryConfigurationParser.cpp
 *
 *  Created on: Apr 8, 2019
 *      Author: jamal
 */

#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/datadriven/datamining/configuration/GeometryConfigurationParser.hpp>

#include <algorithm>
#include <string>

namespace sgpp {
namespace datadriven {

using sgpp::base::data_exception;

StencilType GeometryConfigurationParser::parseStencil(const std::string &input) {
  auto inputLower = input;
  std::transform(inputLower.begin(), inputLower.end(), inputLower.begin(), ::tolower);

  if (inputLower.compare("directneighbour") == 0) {
    return sgpp::datadriven::StencilType::DirectNeighbour;
  } else if(inputLower.compare("allhierarchicalparent") == 0){
    return sgpp::datadriven::StencilType::AllHierarchicalParent;
  } else if(inputLower.compare("nexthierarchicalparent") == 0){
    return sgpp::datadriven::StencilType::NextHierarchicalParent;
  } else if(inputLower.compare("none") == 0){
    return sgpp::datadriven::StencilType::None;
  } else {
    std::string errorMsg = "Failed to convert string \"" + input + "\" to any known "
        "StencilType";
    throw data_exception(errorMsg.c_str());
  }
}
} /* namespace datadriven */
} /* namespace sgpp */