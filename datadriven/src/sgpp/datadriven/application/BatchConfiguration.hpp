// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef BATCHCONFIGURATION_HPP
#define BATCHCONFIGURATION_HPP

#include <sgpp/globaldef.hpp>
#include <string>


//using namespace sg::base;
using namespace std;


namespace SGPP {
  namespace base {
    /**
     * structure to provide parameters for the BatchLearner
     */
    struct BatchConfiguration {
      string filename;//!< arff-file to be read
      size_t batchsize;//!< size of one batch
      size_t samples;//!< number of samles for the monte carlo sampling (normalization) (0=don't sample) good: 1000
      int seed;//!< seed for the sampling
      int wMode;//!< number of weighting mode to use x = batch#, y = wArgument: 0 = all batches are equal, 1 = linear (x*y), 2 = pow(y,x), 3 = y/x, 4 = only the last batch counts, 5 = weigh new batch by proportion, but at least y
      float wArgument;//!< argument for the weighting method
      size_t refineEvery;//!< refine every xth batch (0=never)
      bool verbose;//!< verbose flag
      size_t stack;//!< number of last batches alphavectors to be saved (0=all)
      size_t testsize;//!< how many items to test from the data following the batch (0=don't test after learned)
      float lambda;//!< lambda for solving
    };


  }
}

#endif /* BATCHCONFIGURATION_HPP */
