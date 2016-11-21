// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/optimization/NewtonOptimizer.hpp>

namespace sgpp {
namespace combigrid {

OptimizationGuess NewtonOptimizer::minimize(const OptimizationGuess& guess, size_t maxNumIterations,
                                            double h) {
  OptimizationGuess g = guess;

  for (size_t i = 0; i < maxNumIterations; ++i) {
    double x = g.b - h;
    double y = g.b;
    double z = g.b + h;

    double fx = f(x);
    double fy = g.fb;
    double fz = f(z);

    double secondDeriv = fx + fz - 2 * fy;

    if (secondDeriv <= 0) {  // the function is not locally strictly convex
      return g;
    }

    double firstDeriv = fz - fx;

    if (firstDeriv == 0) {
      return g;
    }

    double xnew = y - firstDeriv * h / (2 * secondDeriv);
    double fxnew = f(xnew);

    if (firstDeriv > 0) {
      // xnew < y
      if (xnew < g.a) {
        xnew = g.a;
      }

      g = OptimizationGuess(g.a, xnew, g.b, g.fa, fxnew, g.fb);
    } else {
      if (xnew > g.c) {
        xnew = g.c;
      }

      g = OptimizationGuess(g.b, xnew, g.c, g.fb, fxnew, g.fc);
    }
  }

  return g;
}

} /* namespace combigrid */
} /* namespace sgpp */
