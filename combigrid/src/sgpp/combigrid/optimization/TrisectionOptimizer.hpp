// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/combigrid/optimization/OptimizationGuess.hpp>

namespace sgpp {
namespace combigrid {

class TrisectionOptimizer {
  SingleFunction f;

 public:
  explicit TrisectionOptimizer(SingleFunction f) : f(f) {}

  OptimizationGuess refine(OptimizationGuess guess, size_t numIterations = 50,
                           size_t numCandidates = 1);
};

} /* namespace combigrid */
} /* namespace sgpp */
