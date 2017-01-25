// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef REGULARIZATIONCONFIGURATION_HPP_
#define REGULARIZATIONCONFIGURATION_HPP_

#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <algorithm>

namespace sgpp {
namespace datadriven {

enum class RegularizationType { Identity, Laplace, Diagonal, Lasso, ElasticNet, GroupLasso };

struct RegularizationConfiguration {
  RegularizationType regType_;
  double lambda_;
  double l1Ratio_;
  double exponentBase_;
};

class RegularizationTypeParser {
 public:
  static RegularizationType parse(const std::string& input) {
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
