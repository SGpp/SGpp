// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/datadriven/configuration/RegularizationConfiguration.hpp>

#include <map>
#include <string>

namespace sgpp {
namespace datadriven {

class RegularizationTypeParser {
 public:
  /**
   * Convert strings to values #sgpp::datadriven::RegularizationType. Throws if there is no valid
   * representation
   * @param input case insensitive string representation of a
   * #sgpp::datadriven::RegularizationType.
   * @return the corresponding #sgpp::datadriven::RegularizationType.
   */
  static RegularizationType parse(const std::string &input);

  /**
   * generate string representations for values of #sgpp::datadriven::RegularizationType.
   * @param type enum value.
   * @return string representation of a #sgpp::datadriven::RegularizationType.
   */
  static const std::string &toString(RegularizationType type);

 private:
  typedef std::map<RegularizationType, std::string> RegularizationTypeMap_t;

  /**
   * Map containing all values of  #sgpp::datadriven::RegularizationType and the corresponding
   * string representation.
   */
  static const RegularizationTypeMap_t regularizationTypeMap;
};
} /* namespace datadriven */
} /* namespace sgpp */
