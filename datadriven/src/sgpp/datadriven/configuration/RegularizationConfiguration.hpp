// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {

enum class RegularizationMetricType { mse, nll, accuracy, residual };

enum class RegularizationType { Identity, Laplace, Diagonal, Lasso, ElasticNet, GroupLasso };

struct RegularizationConfiguration {
  RegularizationType type_ = RegularizationType::Identity;
  double lambda_ = 0.01;
  double l1Ratio_ = 0.0;
  double exponentBase_ = 1.0;
  double lamda_start_ = 0.01;
  double lambda_end_ = 0.01;
  double lambda_steps_ = 0;
  bool lambda_log_scale_ = false;
  bool optimizeLambda_ = false;
  double optimizerTolerance_ = 1e-15;
  double convergenceThreshold_ = 1e-5;
  double intervalA_ = 1e-15;
  double intervalB_ = 1.0;

  RegularizationMetricType regularizationMetric_ = RegularizationMetricType::residual;
};
}  // namespace datadriven
}  // namespace sgpp
