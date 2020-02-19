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
 * Convenience class to convert strings to #sgpp::datadriven::DataSourceFileType and generate
 * string representations for values of #sgpp::datadriven::DataSourceFileType.
 */
class DataSourceFileTypeParser {
 public:
  /**
   * Convert strings to values #sgpp::datadriven::DataSourceFileType. Throws if there is no valid
   * representation
   * @param input case insensitive string representation of a
   * #sgpp::datadriven::DataSourceFileType.
   * @return the corresponding #sgpp::datadriven::DataSourceFileType.
   */
  static DataSourceFileType parse(const std::string &input);

  /**
   * generate string representations for values of #sgpp::datadriven::DataSourceFileType.
   * @param type enum value.
   * @return string representation of a #sgpp::datadriven::DataSourceFileType.
   */
  static const std::string &toString(DataSourceFileType type);

 private:
  using FileTypeMap_t = std::map<DataSourceFileType, std::string>;

  /**
   * Map containing all values of  #sgpp::datadriven::DataSourceFileType and the corresponding
   * string representation.
   */
  static const FileTypeMap_t fileTypeMap;
};
} /* namespace datadriven */
} /* namespace sgpp */
