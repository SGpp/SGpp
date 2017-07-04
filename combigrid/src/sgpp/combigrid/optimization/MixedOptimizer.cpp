// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/optimization/MixedOptimizer.hpp>
#include <sgpp/combigrid/optimization/NewtonOptimizer.hpp>
#include <sgpp/combigrid/optimization/TrisectionOptimizer.hpp>

namespace sgpp {
namespace combigrid {

OptimizationGuess MixedOptimizer::minimize(const OptimizationGuess& guess) {
  TrisectionOptimizer trOpt(f);
  // NewtonOptimizer nOpt(f);

  // double h = (guess.c - guess.a) * 1e-5;

  auto g = trOpt.minimize(guess, 15, 5);
  // g = nOpt.minimize(g, 3, h);

  /*if (g.c - g.b < g.b - g.a) {
    g.a = 2 * g.b - g.c;
    g.fa = f(g.a);
  } else {
    g.c = 2 * g.b - g.a;
    g.fc = f(g.c);
  }*/

  // return trOpt.minimize(guess, 70, 20);
  // return trOpt.minimize(OptimizationGuess::initial(g.b - h, g.b + h, f), 30, 1);
  return trOpt.minimize(g, 20, 1);
}

} /* namespace combigrid */
} /* namespace sgpp */
