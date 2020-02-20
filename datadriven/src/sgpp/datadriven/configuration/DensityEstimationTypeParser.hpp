// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/datadriven/configuration/DensityEstimationConfiguration.hpp>

#include <map>
#include <string>

namespace sgpp {
namespace datadriven {

class DensityEstimationTypeParser {
 public:
  /**
   * Convert strings to values #sgpp::datadriven::DensityEstimationType.
   * Throws if there is no valid representation
   * @param input case insensitive string representation of a
   * #sgpp::datadriven::DensityEstimationType.
   * @return the corresponding #sgpp::datadriven::DensityEstimationType.
   */
  static DensityEstimationType parse(const std::string &input);

  /**
   * generate string representations for values of #sgpp::datadriven::DensityEstimationType.
   * @param type enum value.
   * @return string representation of a #sgpp::datadriven::DensityEstimationType.
   */
  static const std::string &toString(DensityEstimationType type);

 private:
  typedef std::map<DensityEstimationType, std::string> DensityEstimationTypeMap_t;

  /**
   * Map containing all values of  #sgpp::datadriven::DensityEstimationType and the corresponding
   * string representation.
   */
  static const DensityEstimationTypeMap_t densityEstimationTypeMap;
};
} /* namespace datadriven */
} /* namespace sgpp */
