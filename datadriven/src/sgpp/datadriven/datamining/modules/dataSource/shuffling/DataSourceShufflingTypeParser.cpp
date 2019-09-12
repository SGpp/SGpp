// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/shuffling/DataSourceShufflingTypeParser.hpp>

#include <algorithm>
#include <map>
#include <string>

namespace sgpp {
namespace datadriven {

using sgpp::base::data_exception;

DataSourceShufflingType DataSourceShufflingTypeParser::parse(const std::string& input) {
  auto inputLower = input;
  std::transform(inputLower.begin(), inputLower.end(), inputLower.begin(), ::tolower);

  if (inputLower == shufflingTypeMap.at(DataSourceShufflingType::random)) {
    return DataSourceShufflingType::random;
  } else if (inputLower == shufflingTypeMap.at(DataSourceShufflingType::sequential)) {
    return DataSourceShufflingType::sequential;
  } else {
    const auto errorMsg =
        "Failed to convert string \"" + input + "\" to any known ScorerShufflingType";
    throw data_exception(errorMsg.c_str());
  }
}

const std::string& DataSourceShufflingTypeParser::toString(DataSourceShufflingType type) {
  return shufflingTypeMap.at(type);
}

const DataSourceShufflingTypeParser::ShufflingTypeMap_t
  DataSourceShufflingTypeParser::shufflingTypeMap =
    []() {
      return ShufflingTypeMap_t{std::make_pair(DataSourceShufflingType::random, "random"),
                                std::make_pair(DataSourceShufflingType::sequential, "sequential")};
    }();
} /* namespace datadriven */
} /* namespace sgpp */
