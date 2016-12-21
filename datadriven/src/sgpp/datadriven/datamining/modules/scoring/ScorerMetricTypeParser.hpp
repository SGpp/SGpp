/*
 * ScorerMetricTypeParser.hpp
 *
 *  Created on: 21.12.2016
 *      Author: Michael Lettrich
 */

#pragma once

#include <sgpp/datadriven/datamining/modules/scoring/ScorerConfig.hpp>

#include <map>
#include <string>

namespace sgpp {
namespace datadriven {

/**
 * Convenience class to convert strings to #sgpp::datadriven::ScorerMetricType and generate
 * string representations for values of #sgpp::datadriven::ScorerMetricType.
 */
class ScorerMetricTypeParser {
 public:
  /**
   * Convert strings to values #sgpp::datadriven::ScorerMetricType. Throws if there is no valid
   * representation
   * @param input case insensitive string representation of a
   * #sgpp::datadriven::ScorerMetricType.
   * @return the corresponding #sgpp::datadriven::ScorerMetricType.
   */
  static ScorerMetricType parse(const std::string& input);

  /**
   * generate string representations for values of #sgpp::datadriven::ScorerMetricType.
   * @param type enum value.
   * @return string representation of a #sgpp::datadriven::ScorerMetricType.
   */
  static const std::string& stringRepresentation(ScorerMetricType type);

 private:
  /**
   * Map containing all values of  #sgpp::datadriven::ScorerMetricType and the corresponding
   * string representation.
   */
  static const std::map<ScorerMetricType, std::string> metricTypeMap;
};

} /* namespace datadriven */
} /* namespace sgpp */
