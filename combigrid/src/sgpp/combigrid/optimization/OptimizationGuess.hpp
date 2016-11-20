// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/combigrid/SingleFunction.hpp>

namespace sgpp {
namespace combigrid {

class OptimizationGuess {
 public:
  double a, b, c, fa, fb, fc;

  OptimizationGuess(double a, double b, double c, double fa, double fb, double fc)
      : a(a), b(b), c(c), fa(fa), fb(fb), fc(fc) {}

  static OptimizationGuess initial(double a, double c, SingleFunction f) {
    double b = 0.5 * (a + c);
    return OptimizationGuess(a, b, c, f(a), f(b), f(c));
  }
};

bool operator<(OptimizationGuess const &first, OptimizationGuess const &second);

bool operator>(OptimizationGuess const &first, OptimizationGuess const &second);

} /* namespace combigrid */
} /* namespace sgpp */
