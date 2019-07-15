/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * VisualizationTypesParser.hpp
 *
 *  Created on: 09 June, 2019
 *      Author: Vincent Bautista
 */

#pragma once
#include <sgpp/datadriven/datamining/modules/visualization/VisualizationGeneralConfig.hpp>
#include <map>
#include <string>


namespace sgpp {
namespace datadriven {

/**
 * Convenience class to convert strings to #sgpp::datadriven::VisualizationTypesParser and generate
 * string representations for values of #sgpp::datadriven::VisualizationFileTypes and
 * #sgpp::datadriven::VisualizationLibraries
 */
class VisualizationTypesParser {

public:
  /**
   * Convert strings to values #sgpp::datadriven::VisualizationFileTypes. Throws if there is no valid
   * representation
   * @param input case insensitive string representation of a
   * #sgpp::datadriven::VisualizationFileTypes.
   * @return the corresponding #sgpp::datadriven::VisualizationFileTypes.
   */

  static VisualizationFileType parseFileType(const std::string &input);

  /**
   * generate string representations for values of #sgpp::datadriven::DataSourceFileType.
   * @param type enum value.
   * @return string representation of a #sgpp::datadriven::DataSourceFileType.
   */
  static const std::string &toString(VisualizationFileType type);



 private:
  using FileTypeMap_t = std::map<VisualizationFileType, std::string>;

  /**
   * Map containing all values of  #sgpp::datadriven::VisualizationFileType and the corresponding
   * string representation.
   */
  static const FileTypeMap_t fileTypeMap;

};
} /* namespace datadriven */
} /* namespace sgpp */
