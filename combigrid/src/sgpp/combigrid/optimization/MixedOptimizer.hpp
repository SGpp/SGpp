// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/combigrid/optimization/OptimizationGuess.hpp>

namespace sgpp {
namespace combigrid {

class MixedOptimizer {
  SingleFunction f;

 public:
  explicit MixedOptimizer(SingleFunction const &f) : f(f) {}

  OptimizationGuess minimize(OptimizationGuess const &guess);
};

} /* namespace combigrid */
} /* namespace sgpp */
