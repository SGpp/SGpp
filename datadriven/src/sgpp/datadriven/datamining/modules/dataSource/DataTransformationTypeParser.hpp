// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/datadriven/datamining/modules/dataSource/DataTransformationConfig.hpp>

#include <map>
#include <string>

namespace sgpp {
namespace datadriven {

/**
 * Convenience class to convert strings to #sgpp::datadriven::DataTransformationType and generate
 * string representations for values of #sgpp::datadriven::DataTransformationType.
 */
class DataTransformationTypeParser {
 public:
  /**
   * Convert strings to values #sgpp::datadriven::DataTransformationType. Throws if there is no valid
   * representation
   * @param input case insensitive string representation of a #sgpp::datadriven::DataTransformationType.
   * @return the corresponding #sgpp::datadriven::DataTransformationType.
   */
  static DataTransformationType parse(const std::string &input);

  /**
   * Generate string representations for values of #sgpp::datadriven::DataTransformationType.
   * @param type enum value.
   * @return string representation of a #sgpp::datadriven::DataTransformationType.
   */
  static const std::string &toString(DataTransformationType type);

 private:
  typedef std::map<DataTransformationType, std::string> TransformationTypeMap_t;

  /**
   * Map containing all values of  #sgpp::datadriven::DataTransformationType and the corresponding
   * string representation.
   */
  static const TransformationTypeMap_t transformationTypeMap;
};
} /* namespace datadriven */
} /* namespace sgpp */
