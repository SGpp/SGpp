/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * RegularizationTypeParser.hpp
 *
 * Created on: Jan 25, 2017
 *     Author: Michael Lettrich
 */

#pragma once

#include <sgpp/datadriven/application/RegularizationConfiguration.hpp>

#include <map>
#include <string>

namespace sgpp {
namespace datadriven {

class RegularizationTypeParser {
 public:
  /**
   * Convert strings to values #sgpp::datadriven::DataSourceFileType. Throws if there is no valid
   * representation
   * @param input case insensitive string representation of a
   * #sgpp::datadriven::DataSourceFileType.
   * @return the corresponding #sgpp::datadriven::DataSourceFileType.
   */
  static RegularizationType parse(const std::string& input);

  /**
   * generate string representations for values of #sgpp::datadriven::DataSourceFileType.
   * @param type enum value.
   * @return string representation of a #sgpp::datadriven::DataSourceFileType.
   */
  static const std::string& toString(RegularizationType type);

 private:
  using RegularizationTypeMap_t = std::map<RegularizationType, std::string>;

  /**
   * Map containing all values of  #sgpp::datadriven::DataSourceFileType and the corresponding
   * string representation.
   */
  static const RegularizationTypeMap_t regularizationTypeMap;
};

} /* namespace datadriven */
} /* namespace sgpp */
