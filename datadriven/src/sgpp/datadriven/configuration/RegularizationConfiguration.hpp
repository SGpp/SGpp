// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef REGULARIZATIONCONFIGURATION_HPP_
#define REGULARIZATIONCONFIGURATION_HPP_

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {

enum class RegularizationType { Identity, Laplace, Diagonal, Lasso, ElasticNet, GroupLasso };

struct RegularizationConfiguration {
  RegularizationType type_ = RegularizationType::Identity;
  double lambda_ = 0.01;
  double l1Ratio_ = 0.0;
  double exponentBase_ = 1.0;
};
}  // namespace datadriven
}  // namespace sgpp

#endif /* REGULARIZATIONCONFIGURATION_HPP_ */
