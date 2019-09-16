// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/datadriven/datamining/modules/dataSource/DataSourceConfig.hpp>

#include <map>
#include <string>

namespace sgpp {
namespace datadriven {

/**
 * Convenience class to convert strings to #sgpp::datadriven::DataSourceShufflingType and generate
 * string representations for values of #sgpp::datadriven::DataSourceShufflingType.
 */
class DataSourceShufflingTypeParser {
 public:
  /**
   * Convert strings to values #sgpp::datadriven::DataSourceShufflingType. Throws if there is no valid
   * representation
   * @param input case insensitive string representation of a
   * #sgpp::datadriven::DataSourceShufflingType.
   * @return the corresponding #sgpp::datadriven::DataSourceShufflingType.
   */
  static DataSourceShufflingType parse(const std::string& input);

  /**
   * generate string representations for values of #sgpp::datadriven::DataSourceShufflingType.
   * @param type enum value.
   * @return string representation of a #sgpp::datadriven::DataSourceShufflingType.
   */
  static const std::string& toString(DataSourceShufflingType type);

 private:
  typedef std::map<DataSourceShufflingType, std::string> ShufflingTypeMap_t;
  /**
   * Map containing all values of  #sgpp::datadriven::DataSourceShufflingType and the corresponding
   * string representation.
   */
  static const ShufflingTypeMap_t shufflingTypeMap;
};

} /* namespace datadriven */
} /* namespace sgpp */
