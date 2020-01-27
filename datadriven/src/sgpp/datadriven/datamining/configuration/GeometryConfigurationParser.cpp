// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

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
  } else if (inputLower.compare("allhierarchicalparent") == 0) {
    return sgpp::datadriven::StencilType::AllHierarchicalParent;
  } else if (inputLower.compare("nexthierarchicalparent") == 0) {
    return sgpp::datadriven::StencilType::NextHierarchicalParent;
  } else if (inputLower.compare("block") == 0) {
    return sgpp::datadriven::StencilType::Block;
  } else if (inputLower.compare("none") == 0) {
    return sgpp::datadriven::StencilType::None;
  } else {
    std::string errorMsg = "Failed to convert string \"" + input +
                           "\" to any known "
                           "StencilType";
    throw data_exception(errorMsg.c_str());
  }
}

const std::string &GeometryConfigurationParser::toString(StencilType type) { return stencilTypeMap.at(type); }

const GeometryConfigurationParser::StencilTypeMap_t GeometryConfigurationParser::stencilTypeMap = []() {
  return StencilTypeMap_t{
      std::make_pair(StencilType::DirectNeighbour, "DirectNeighbour"),
      std::make_pair(StencilType::AllHierarchicalParent, "AllHierarchicalParent"),
      std::make_pair(StencilType::NextHierarchicalParent, "NextHierarchicalParent"),
      std::make_pair(StencilType::Block, "Block"),
      std::make_pair(StencilType::None, "None")};
}();
} /* namespace datadriven */
} /* namespace sgpp */
