// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/combigrid/optimization/OptimizationGuess.hpp>

namespace sgpp {
namespace combigrid {

class NewtonOptimizer {
  SingleFunction f;

 public:
  explicit NewtonOptimizer(SingleFunction const &f) : f(f) {}

  OptimizationGuess minimize(OptimizationGuess const &guess, size_t maxNumIterations,
                             double h = 1e-6);
};

} /* namespace combigrid */
} /* namespace sgpp */
