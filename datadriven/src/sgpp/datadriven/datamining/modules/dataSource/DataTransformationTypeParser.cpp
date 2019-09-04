// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/modules/dataSource/DataTransformationTypeParser.hpp>

#include <algorithm>
#include <string>

namespace sgpp {
namespace datadriven {

DataTransformationType DataTransformationTypeParser::parse(const std::string &input) {
  auto inputLower = input;
  std::transform(inputLower.begin(), inputLower.end(), inputLower.begin(), ::tolower);

  if (inputLower.compare("rosenblatt") == 0) {
    return DataTransformationType::ROSENBLATT;
  } else {
    return DataTransformationType::NONE;
  }
}

const std::string &DataTransformationTypeParser::toString(DataTransformationType type) {
  return transformationTypeMap.at(type);
}

const DataTransformationTypeParser::TransformationTypeMap_t
    DataTransformationTypeParser::transformationTypeMap = []() {
  return DataTransformationTypeParser::TransformationTypeMap_t{
      std::make_pair(DataTransformationType::NONE, "None"),
      std::make_pair(DataTransformationType::ROSENBLATT, "Rosenblatt")};
}();
} /* namespace datadriven */
} /* namespace sgpp */
