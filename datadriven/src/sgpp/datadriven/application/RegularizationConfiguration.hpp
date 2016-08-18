// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef REGULARIZATIONCONFIGURATION_HPP_
#define REGULARIZATIONCONFIGURATION_HPP_

#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/globaldef.hpp>

#include <algorithm>

namespace sgpp {
namespace datadriven {

enum class RegularizationType { Identity, Laplace };

struct RegularizationConfiguration {
  RegularizationType regType_;
};

class RegularizationTypeParser {
 public:
  RegularizationType operator()(const std::string& input) const {
    auto inputLower = input;
    std::transform(inputLower.begin(), inputLower.end(), inputLower.begin(), ::tolower);

    if (inputLower.compare("identity") == 0) {
      return sgpp::datadriven::RegularizationType::Identity;
    } else if (inputLower.compare("laplace") == 0) {
      return sgpp::datadriven::RegularizationType::Laplace;
    } else {
      std::string errorMsg =
          "Failed to convert string \"" + input + "\" to any known RegularizationType";
      throw base::data_exception(errorMsg.c_str());
    }
  }
};

}  // namespace datadriven
}  // namespace sgpp

#endif /* REGULARIZATIONCONFIGURATION_HPP_ */
