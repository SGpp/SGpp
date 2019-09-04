// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/datadriven/datamining/modules/fitting/FitterConfiguration.hpp>

#include <map>
#include <string>

namespace sgpp {
namespace datadriven {

/**
 * Convenience class to convert strings to #sgpp::datadriven::FitterType and generate
 * string representations for values of #sgpp::datadriven::FitterType.
 */
class FitterTypeParser {
 public:
  /**
   * Convert strings to values #sgpp::datadriven::FitterType. Throws if there is no valid
   * representation
   * @param input case insensitive string representation of a
   * #sgpp::datadriven::FitterType.
   * @return the corresponding #sgpp::datadriven::FitterType.
   */
  static FitterType parse(const std::string &input);

  /**
   * generate string representations for values of #sgpp::datadriven::ScorerMetricType.
   * @param type enum value.
   * @return string representation of a #sgpp::datadriven::FitterType.
   */
  static const std::string &toString(FitterType type);

 private:
  typedef std::map<FitterType, std::string> FitterTypeMap_t;

  /**
   * Map containing all values of #sgpp::datadriven::FitterType and the corresponding
   * string representation.
   */
  static const FitterTypeMap_t fitterTypeMap;
};
} /* namespace datadriven */
} /* namespace sgpp */
