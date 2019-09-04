// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>

/**
 * Struct that stores all the configuration information for crossvalidation
 */

namespace sgpp {
namespace datadriven {

struct CrossvalidationConfiguration {
  // parameters for cross-validation
  bool enable_ = false;   // enables cross-validation
  size_t kfold_ = 5;  // number of batches for cross validation
  int seed_;      // seed for randomized k-fold
  bool shuffle_;  // randomized/sequential k-fold
  bool silent_;   // verbosity

  // regularization parameter optimization
  double lambda_;       // regularization parameter
  double lambdaStart_;  // lower bound for lambda search range
  double lambdaEnd_;    // upper bound for lambda search range
  // number of lambdas to be tested within the range defined by lambdaStart and
  // lambdaEdns;
  // must be > 1
  size_t lambdaSteps_;
  bool logScale_;  // search the optimization interval on a log-scale
};

}  // namespace datadriven
}  // namespace sgpp
