// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

namespace sgpp {
namespace datadriven {

/**
 * Structure that contains information about the learners behaviour
 */
struct LearnerConfiguration {
  /**
  * Weigting factor for older batches
  */
  double learningRate_ = 1.0;

  /**
   * Determine if the relative frequency of instances of a class should be used as prior
   * (false corresponds to a uniform prior)
   */
  bool usePrior_ = false;
};
}  // namespace datadriven
}  // namespace sgpp
