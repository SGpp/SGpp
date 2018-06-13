/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * ScorerConfig.hpp
 *
 * Created on: Aug 18, 2016
 *     Author: Michael Lettrich
 */

#pragma once

#include <string>

namespace sgpp {
namespace datadriven {
/**
 * Enumeration of all supported shuffling types used to permute samples in a dataset. An entry
 * exists for each object that derives from #sgpp::datadriven::ShufflingFunctor. Used for
 * configuration and factory methods.
 */
enum class ScorerShufflingType { random, sequential };
/**
 * Enumeration of all supported metrics used to quantify approximation quality of a trained model.
 * An entry exists for each object that derives from #sgpp::datadriven::Metric. Used for
 * configuration and factory methods.
 */
enum class ScorerMetricType { mse, nll };

/**
 * Set of parameters required to fully configure #sgpp::datadriven::CrossValidation objects.
 */
struct CrossValidationConfiguration {
  /**
   * Number of folds used for cross validation.
   */
  size_t folds = 5;
  /**
   * Random seed used for permutation of samples in the shuffling functor. For each random seed,
   * the same set of random numbers are generated, meaning that same random seeds generate the same
   * sequence of random numbers. Default value of -1 indicates to generate an arbitrary random seed.
   */
  int64_t randomSeed = -1;
  /**
   * Type of metric that should be used to calculate the accuracy of the fit.
   */
  ScorerMetricType metric;
  /**
   * Type of shuffling operator to be used for rearrangement of samples in the parsed dataset.
   */
  ScorerShufflingType shuffling = ScorerShufflingType::random;
};

/**
 * Set of parameters required to fully configure #sgpp::datadriven::SplittingScorer objects.
 */
struct TestingConfiguration {
  /**
   * value between 0 and 1 to specify the ratio between testing set and training set.
   */
  double testingPortion = 0.3;
  /**
   * Random seed used for permutation of samples in the shuffling functor. For each random seed,
   * the same set of random numbers are generated, meaning that same random seeds generate the same
   * sequence of random numbers. Default value of -1 indicates to generate an arbitrary random seed.
   */
  int64_t randomSeed = -1;
  /**
   * Type of metric that should be used to calculate the accuracy of the fit.
   */
  ScorerMetricType metric;
  /**
   * Type of shuffling operator to be used for rearrangement of samples in the parsed dataset.
   */
  ScorerShufflingType shuffling = ScorerShufflingType::random;
};

} /* namespace datadriven */
} /* namespace sgpp */
