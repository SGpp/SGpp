/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * ScorerShufflingTypeParser.cpp
 *
 *  Created on: 21.12.2016
 *      Author: Michael Lettrich
 */

#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/ScorerShufflingTypeParser.hpp>

#include <algorithm>
#include <map>
#include <string>

namespace sgpp {
namespace datadriven {

using sgpp::base::data_exception;

ScorerShufflingType ScorerShufflingTypeParser::parse(const std::string &input) {
  auto inputLower = input;
  std::transform(inputLower.begin(), inputLower.end(), inputLower.begin(), ::tolower);

  if (inputLower == shufflingTypeMap.at(ScorerShufflingType::random)) {
    return ScorerShufflingType::random;
  } else if (inputLower == shufflingTypeMap.at(ScorerShufflingType::sequential)) {
    return ScorerShufflingType::sequential;
  } else {
    const auto errorMsg =
        "Failed to convert string \"" + input + "\" to any known ScorerShufflingType";
    throw data_exception(errorMsg.c_str());
  }
}

const std::string &ScorerShufflingTypeParser::toString(ScorerShufflingType type) {
  return shufflingTypeMap.at(type);
}

const ScorerShufflingTypeParser::ShufflingTypeMap_t ScorerShufflingTypeParser::shufflingTypeMap =
    []() {
      return ShufflingTypeMap_t{std::make_pair(ScorerShufflingType::random, "random"),
                                std::make_pair(ScorerShufflingType::sequential, "sequential")};
    }();
} /* namespace datadriven */
} /* namespace sgpp */
