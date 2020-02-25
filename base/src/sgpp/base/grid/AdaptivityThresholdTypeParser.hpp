// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/grid/Grid.hpp>

#include <map>
#include <string>

namespace sgpp {
namespace base {
using sgpp::base::AdaptivityThresholdType;

class AdaptivityThresholdTypeParser {
 public:
  /**
   * Convert strings to values #sgpp::base::AdaptivityThresholdType. Throws if there is no valid
   * representation
   * @param input case insensitive string representation of a
   * #sgpp::base::AdaptivityThresholdType.
   * @return the corresponding #sgpp::base::AdaptivityThresholdType.
   */
  static AdaptivityThresholdType parse(const std::string& input);

  /**
   * generate string representations for values of #sgpp::base::AdaptivityThresholdType.
   * @param type enum value.
   * @return string representation of a #sgpp::base::AdaptivityThresholdType.
   */
  static const std::string& toString(AdaptivityThresholdType type);

 private:
  typedef std::map<AdaptivityThresholdType, std::string> AdaptivityThresholdTypeMap_t;

  /**
   * Map containing all values of #sgpp::base::AdaptivityThresholdType and the corresponding
   * string representation.
   */
  static const AdaptivityThresholdTypeMap_t refinementFunctorTypeMap;
};

} /* namespace base */
} /* namespace sgpp */
