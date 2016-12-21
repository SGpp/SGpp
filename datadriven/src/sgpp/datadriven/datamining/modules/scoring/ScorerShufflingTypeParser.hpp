/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * ScorerShufflingTypeParser.hpp
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
 * Convenience class to convert strings to #sgpp::datadriven::ScorerShufflingType and generate
 * string representations for values of #sgpp::datadriven::ScorerShufflingType.
 */
class ScorerShufflingTypeParser {
 public:
  /**
   * Convert strings to values #sgpp::datadriven::ScorerShufflingType. Throws if there is no valid
   * representation
   * @param input case insensitive string representation of a
   * #sgpp::datadriven::ScorerShufflingType.
   * @return the corresponding #sgpp::datadriven::ScorerShufflingType.
   */
  static ScorerShufflingType parse(const std::string& input);

  /**
   * generate string representations for values of #sgpp::datadriven::ScorerShufflingType.
   * @param type enum value.
   * @return string representation of a #sgpp::datadriven::ScorerShufflingType.
   */
  static const std::string& stringRepresentation(ScorerShufflingType type);

 private:
  /**
   * Map containing all values of  #sgpp::datadriven::ScorerShufflingType and the corresponding
   * string representation.
   */
  static const std::map<ScorerShufflingType, std::string> shufflingTypeMap;
};

} /* namespace datadriven */
} /* namespace sgpp */
