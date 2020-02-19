// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/GeneralGridTypeParser.hpp>

#include <sgpp/base/exception/data_exception.hpp>
#include <algorithm>
#include <string>

namespace sgpp {
namespace base {

using sgpp::base::data_exception;

sgpp::base::GeneralGridType GeneralGridTypeParser::parse(
    const std::string &input) {
  auto inputLower = input;
  std::transform(inputLower.begin(), inputLower.end(), inputLower.begin(),
                 ::tolower);
  if (inputLower == "regular") {
    return sgpp::base::GeneralGridType::RegularSparseGrid;
  } else if (inputLower == "refinedcoarsened") {
    return sgpp::base::GeneralGridType::RefinedCoarsenedSparseGrid;
  } else if (inputLower == "geometryaware") {
    return sgpp::base::GeneralGridType::GeometryAwareSparseGrid;
  } else if (inputLower == "component") {
    return sgpp::base::GeneralGridType::ComponentGrid;
  } else {
    std::string what = "Failed to parse general grid type \"" + input + "\".";
    throw data_exception(what.c_str());
  }
}

const std::string &GeneralGridTypeParser::toString(
    sgpp::base::GeneralGridType type) {
  return generalGridTypeMap.at(type);
}

const GeneralGridTypeParser::GeneralGridTypeMap_t
    GeneralGridTypeParser::generalGridTypeMap = []() {
      return GeneralGridTypeMap_t{
          std::make_pair(sgpp::base::GeneralGridType::RegularSparseGrid,
                         "regular"),
          std::make_pair(
              sgpp::base::GeneralGridType::RefinedCoarsenedSparseGrid,
              "refinedcoarsened"),
          std::make_pair(sgpp::base::GeneralGridType::GeometryAwareSparseGrid,
                         "geometryaware"),
          std::make_pair(sgpp::base::GeneralGridType::ComponentGrid,
                         "component")};
    }();
} /* namespace base */
} /* namespace sgpp */
