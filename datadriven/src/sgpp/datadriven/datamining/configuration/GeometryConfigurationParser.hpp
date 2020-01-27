// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/datadriven/configuration/GeometryConfiguration.hpp>

#include <map>
#include <string>

namespace sgpp {
namespace datadriven {

class GeometryConfigurationParser {
 public:
  /**
   * Convert strings to values #sgpp::datadriven::StencilType.
   * Throws if there is no valid representation
   * @param input case insensitive string representation of a
   * #sgpp::datadriven::StencilType.
   * @return the corresponding #sgpp::datadriven::StencilType.
   */
  static StencilType parseStencil(const std::string &input);

  /**
   * generate string representations for values of #sgpp::base::StencilType.
   * @param type enum value.
   * @return string representation of a #sgpp::base::StencilType.
   */
  static const std::string &toString(StencilType type);

 private:
  typedef std::map<StencilType, std::string> StencilTypeMap_t;

  /**
   * Map containing all values of  #sgpp::base::StencilType and the corresponding
   * string representation.
   */
  static const StencilTypeMap_t stencilTypeMap;
};
} /* namespace datadriven */
} /* namespace sgpp */
