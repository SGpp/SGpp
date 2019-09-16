// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/modules/dataSource/DataSourceFileTypeParser.hpp>

#include <sgpp/base/exception/data_exception.hpp>

#include <algorithm>
#include <string>

namespace sgpp {
namespace datadriven {

using sgpp::base::data_exception;

DataSourceFileType DataSourceFileTypeParser::parse(const std::string &input) {
  auto inputLower = input;
  std::transform(inputLower.begin(), inputLower.end(), inputLower.begin(), ::tolower);

  if (inputLower == "arff") {
    return DataSourceFileType::ARFF;
  } else if (inputLower == "none") {
    return DataSourceFileType::NONE;
  } else if (inputLower == "csv") {
    return DataSourceFileType::CSV;
  } else {
    const std::string errorMsg =
        "Failed to convert string \"" + input + "\" to any known DataSourceFileType";
    throw data_exception(errorMsg.c_str());
  }
}

const std::string &sgpp::datadriven::DataSourceFileTypeParser::toString(DataSourceFileType type) {
  return fileTypeMap.at(type);
}

const DataSourceFileTypeParser::FileTypeMap_t DataSourceFileTypeParser::fileTypeMap = []() {
  return DataSourceFileTypeParser::FileTypeMap_t{std::make_pair(DataSourceFileType::NONE, "None"),
                                                 std::make_pair(DataSourceFileType::ARFF, "ARFF"),
                                                 std::make_pair(DataSourceFileType::CSV, "CSV")};
}();
} /* namespace datadriven */
} /* namespace sgpp */
