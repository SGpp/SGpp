// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <cstdint>
#include <vector>

namespace sgpp {
namespace datadriven {

/**
 * Configuration for HyperparameterOptimizer
 */
class HPOConfig {
 public:
  HPOConfig() = default;

  void setupDefaults();

  int64_t getSeed() const;

  void setSeed(int64_t seed);

  const std::vector<int64_t> &getStages() const;

  void setStages(const std::vector<int64_t> &stages);

  const std::vector<int64_t> &getConstraints() const;

  void setConstraints(const std::vector<int64_t> &constraints);

  double getLambda() const;

  void setLambda(double lambda);

  int64_t getNRandom() const;

  void setNRandom(int64_t nRandom);

  int64_t getNRuns() const;

  void setNRuns(int64_t nRuns);

  int64_t getNTrainSamples() const;

  void setNTrainSamples(int64_t nTrainSamples);

 private:
  /**
   * Seed for random sampling in both harmonica and bayesian optimization
   */
  int64_t seed;
  /**
   * Number of Training Samples used for HPO
   */
  int64_t nTrainSamples;
  /**
   * Amount of samples to take in each stage of harmonica
   */
  std::vector<int64_t> stages;
  /**
   * Amount of constraints to introduce after each stage of harmonica
   */
  std::vector<int64_t> constraints;
  /**
   * Regularization Lambda used for Lasso Regression in harmonica
   */
  double lambda;
  /**
   * Number of Random samples used to warm up bayesian optimization
   */
  int64_t nRandom;
  /**
   * number of samples bayesian optimization is run for
   */
  int64_t nRuns;
};
} /* namespace datadriven */
} /* namespace sgpp */
