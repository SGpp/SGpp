// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/modules/visualization/VisualizationTypesParser.hpp>

#include <sgpp/base/exception/data_exception.hpp>

#include <algorithm>
#include <string>

namespace sgpp {
namespace datadriven {

using sgpp::base::data_exception;

VisualizationFileType VisualizationTypesParser::parseFileType(const std::string &input) {
  auto inputLower = input;
  std::transform(inputLower.begin(), inputLower.end(), inputLower.begin(), ::tolower);

  if (inputLower == "json") {
    return VisualizationFileType::json;
  } else if (inputLower == "csv") {
    return VisualizationFileType::CSV;
  } else {
    const std::string errorMsg =
      "Failed to convert string \"" + input + "\" to any known VisualizationFileType";
    throw data_exception(errorMsg.c_str());
}
}

const std::string &sgpp::datadriven::VisualizationTypesParser::toString
(VisualizationFileType type) {
  return fileTypeMap.at(type);
}


const VisualizationTypesParser::FileTypeMap_t VisualizationTypesParser::fileTypeMap = []() {
  return VisualizationTypesParser::FileTypeMap_t {
    std::make_pair(VisualizationFileType::json, "ARFF"),
    std::make_pair(VisualizationFileType::CSV, "CSV") };
}();


} /* namespace datadriven */
} /* namespace sgpp */
