// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <string>

namespace sgpp {
namespace datadriven {

/**
 * Configuration structure for Rosenblatt transformation including default values.
 */
struct RosenblattTransformationConfig {
  /**
   * Parameters for calculation of PDF / density estimation
   */
  size_t numSamples_ = 1000;

  size_t gridLevel_ = 2;

  size_t solverMaxIterations_ = 1000;
  double solverEps_ = 1e-10;
  double solverThreshold_ = 1e-10;
};
} /* namespace datadriven */
} /* namespace sgpp */
