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
  RegularizationType type_;
  double lambda_;
  double lamda_start_;
  double lambda_end_;
  double lambda_steps_;
  bool lambda_log_scale_;
  double l1Ratio_;
  double exponentBase_;
  bool optimizeLambda_;
  double optimizerTolerance_;
  double convergenceThreshold_;
  double intervalA_;
  double intervalB_;
};
}  // namespace datadriven
}  // namespace sgpp

#endif /* REGULARIZATIONCONFIGURATION_HPP_ */
