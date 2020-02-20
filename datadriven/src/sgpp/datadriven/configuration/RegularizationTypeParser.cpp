// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/configuration/RegularizationTypeParser.hpp>
#include <sgpp/base/exception/data_exception.hpp>

#include <algorithm>
#include <string>

namespace sgpp {
namespace datadriven {
RegularizationType RegularizationTypeParser::parse(const std::string &input) {
  auto inputLower = input;
  std::transform(inputLower.begin(), inputLower.end(), inputLower.begin(), ::tolower);

  if (inputLower.compare("identity") == 0) {
    return sgpp::datadriven::RegularizationType::Identity;
  } else if (inputLower.compare("laplace") == 0) {
    return sgpp::datadriven::RegularizationType::Laplace;
  } else if (inputLower.compare("diagonal") == 0) {
    return sgpp::datadriven::RegularizationType::Diagonal;
  } else if (inputLower.compare("lasso") == 0) {
    return sgpp::datadriven::RegularizationType::Lasso;
  } else if (inputLower.compare("elasticnet") == 0) {
    return sgpp::datadriven::RegularizationType::ElasticNet;
  } else if (inputLower.compare("grouplasso") == 0) {
    return sgpp::datadriven::RegularizationType::GroupLasso;
  } else {
    std::string errorMsg =
        "Failed to convert string \"" + input + "\" to any known RegularizationType";
    throw base::data_exception(errorMsg.c_str());
  }
}
const std::string &RegularizationTypeParser::toString(RegularizationType type) {
  return regularizationTypeMap.at(type);
}

const RegularizationTypeParser::RegularizationTypeMap_t
    RegularizationTypeParser::regularizationTypeMap = []() {
  return RegularizationTypeParser::RegularizationTypeMap_t{
      std::make_pair(RegularizationType::Identity, "Identity"),
      std::make_pair(RegularizationType::Laplace, "Laplace"),
      std::make_pair(RegularizationType::Diagonal, "Diagonal"),
      std::make_pair(RegularizationType::Lasso, "Lasso"),
      std::make_pair(RegularizationType::ElasticNet, "ElasticNet"),
      std::make_pair(RegularizationType::GroupLasso, "GroupLasso")};
}();
} /* namespace datadriven */
} /* namespace sgpp */
