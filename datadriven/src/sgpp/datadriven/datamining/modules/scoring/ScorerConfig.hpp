// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <string>

namespace sgpp {
namespace datadriven {
/**
 * Enumeration of all supported metrics used to quantify approximation quality of a trained model.
 * An entry exists for each object that derives from #sgpp::datadriven::Metric. Used for
 * configuration and factory methods.
 */
enum class ScorerMetricType { mse, nll, accuracy };

/**
 * Set of parameters to define a scorer instance
 */
struct ScorerConfiguration {
  /**
   * Type of metric that should be used to calculate the accuracy of the fit.
   */
  ScorerMetricType metric_ = ScorerMetricType::accuracy;
};
} /* namespace datadriven */
} /* namespace sgpp */
