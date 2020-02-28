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
  size_t numSamples = 1000;

  size_t gridLevel = 2;

  size_t solverMaxIterations = 1000;
  double solverEps = 1e-10;
  double solverThreshold = 1e-10;
};
} /* namespace datadriven */
} /* namespace sgpp */
