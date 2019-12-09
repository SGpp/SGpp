// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <string>
#include <vector>


namespace sgpp {
namespace datadriven {

/**
 * Configuration structure for Normalization transformation including default values.
 */
struct NormalizationTransformationConfig {
  /**
   * Parameters for setting up the transformation with default values
   */
  std::string method = "minmax";
  bool manualInput = false;
  std::vector<std::vector<double>>  minmaxInput = {};
  unsigned int searchInstances = 100;
  bool stdDeviationHeuristic = true;
  double deviationHeuristic = 3.5;
  unsigned int minmaxStdDeviation = 3;
  std::string outOfBound = "abort";
};
} /* namespace datadriven */
} /* namespace sgpp */
