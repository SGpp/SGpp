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
  size_t kfold_ = 5;      // number of batches for cross validation
  int seed_ = 0;          // seed for randomized k-fold
  bool shuffle_ = false;  // randomized/sequential k-fold
  bool silent_ = true;    // verbosity

  // regularization parameter optimization
  double lambda_ = 1e-3;       // regularization parameter
  double lambdaStart_ = 1e-3;  // lower bound for lambda search range
  double lambdaEnd_ = 1e-3;    // upper bound for lambda search range
  // number of lambdas to be tested within the range defined by lambdaStart and
  // lambdaEdns;
  // must be > 1
  size_t lambdaSteps_ = 0;
  bool logScale_ = false;  // search the optimization interval on a log-scale

  /*
  // Debug method to neatly print internal data
  void dumpToStream(std::ostream& stream_out = std::cout) const {
    stream_out << "enable \t\t\t" << std::boolalpha << enable_ << std::endl;
    stream_out << "kfold \t\t\t" << kfold_ << std::endl;
    stream_out << "seed \t\t\t" << seed_ << std::endl;
    stream_out << "shuffle \t\t" << std::boolalpha << shuffle_ << std::endl;
    stream_out << "silent \t\t\t" << std::boolalpha << silent_ << std::endl;
    stream_out << "lambda \t\t\t" << lambda_ << std::endl;
    stream_out << "lambdaEnd \t\t" << lambdaEnd_ << std::endl;
    stream_out << "lambdaStart \t\t" << lambdaStart_ << std::endl;
    stream_out << "lambdaSteps \t\t" << lambdaSteps_ << std::endl;
    stream_out << "logScale \t\t" << std::boolalpha << logScale_ << std::endl;
  }
  */
};

}  // namespace datadriven
}  // namespace sgpp
