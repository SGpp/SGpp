// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef BATCHCONFIGURATION_HPP
#define BATCHCONFIGURATION_HPP

#include <sgpp/globaldef.hpp>
#include <string>


namespace sgpp {
namespace base {
/**
 * structure to provide parameters for the BatchLearner
 */
struct BatchConfiguration {
  std::string filename_;  //!< arff-file to be read
  size_t batchsize_;  //!< size of one batch
  //!< number of samles for the monte carlo sampling (normalization) (0=don't sample) good: 1000
  size_t samples_;
  int seed_;  //!< seed for the sampling
  //!< number of weighting mode to use x = batch#, y = wArgument: 0 = all batches are equal,
  // 1 = linear (x*y), 2 = pow(y,x), 3 = y/x, 4 = only the last batch counts, 5 = weight new
  // batch by proportion, but at least y
  int wMode_;
  double wArgument_;  //!< argument for the weighting method
  size_t refineEvery_;  //!< refine every xth batch (0=never)
  bool verbose_;  //!< verbose flag
  size_t stack_;  //!< number of last batches alphavectors to be saved (0=all)
  //!< how many items to test from the data following the batch (0=don't test after learned)
  size_t testsize_;
  double lambda_;  //!< lambda for solving
};


}  // namespace base
}  // namespace sgpp

#endif /* BATCHCONFIGURATION_HPP */
