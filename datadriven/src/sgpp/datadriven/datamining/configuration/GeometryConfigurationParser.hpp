/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * GeometryConfigurationParser.hpp
 *
 *  Created on: Apr 8, 2019
 *      Author: jamal
 */

#pragma once

#include <sgpp/datadriven/configuration/GeometryConfiguration.hpp>

#include <map>
#include <string>

namespace sgpp {
namespace datadriven {

class GeometryConfigurationParser {
 public:
  /**
   * Convert strings to values #sgpp::datadriven::DensityEstimationType.
   * Throws if there is no valid representation
   * @param input case insensitive string representation of a
   * #sgpp::datadriven::DensityEstimationType.
   * @return the corresponding #sgpp::datadriven::DensityEstimationType.
   */
  static StencilType parse(const std::string &input);
};
} /* namespace datadriven */
} /* namespace sgpp */
