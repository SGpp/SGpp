// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/optimization/TrisectionOptimizer.hpp>
#include <sgpp/combigrid/utils/BinaryHeap.hpp>

#include <algorithm>
#include <functional>
#include <memory>

namespace sgpp {
namespace combigrid {

OptimizationGuess TrisectionOptimizer::minimize(OptimizationGuess guess, size_t numIterations,
                                              size_t numCandidates) {
  auto oldHeap = std::make_shared<BinaryHeap<OptimizationGuess, std::greater<OptimizationGuess>>>();
  auto newHeap = std::make_shared<BinaryHeap<OptimizationGuess, std::greater<OptimizationGuess>>>();

  oldHeap->push(guess);

  for (size_t i = 0; i < numIterations; ++i) {
    for (size_t j = 0; j < numCandidates && !oldHeap->empty(); ++j) {
      auto g = oldHeap->top();
      oldHeap->pop();
      double ab = 0.5 * (g.a + g.b);
      double bc = 0.5 * (g.b + g.c);
      double fab = f(ab);
      double fbc = f(bc);
      newHeap->push(OptimizationGuess(g.a, ab, g.b, g.fa, fab, g.fb));
      newHeap->push(OptimizationGuess(ab, g.b, bc, fab, g.fb, fbc));
      newHeap->push(OptimizationGuess(g.b, bc, g.c, g.fb, fbc, g.fc));
    }

    std::swap(oldHeap, newHeap);
    newHeap->clear();
  }

  return oldHeap->top();
}

} /* namespace combigrid */
} /* namespace sgpp */
