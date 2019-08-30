// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/operation/multidim/AdaptiveRefinementStrategy.hpp>

#include <stdexcept>
#include <cmath>
#include <vector>

namespace sgpp {
namespace combigrid {

AdaptiveRefinementStrategy::AdaptiveRefinementStrategy(priority_function func) : func(func) {}

double AdaptiveRefinementStrategy::computePriority(std::vector<double> const &predecessorNorms,
                                                   size_t numNewPoints) {
  return func(predecessorNorms, numNewPoints);
}

AdaptiveRefinementStrategy AdaptiveRefinementStrategy::maxStrategy() {
  return AdaptiveRefinementStrategy(
      [](std::vector<double> const &predecessorNorms, size_t numNewPoints) -> double {
        if (predecessorNorms.size() == 0)
          throw std::runtime_error("No predecessors in AdaptiveRefinementStrategy");

        double val = predecessorNorms[0];

        for (size_t i = 1; i < predecessorNorms.size(); ++i) {
          double other = predecessorNorms[i];
          if (other > val) {
            val = other;
          }
        }

        return val / static_cast<double>(numNewPoints);
      });
}

AdaptiveRefinementStrategy AdaptiveRefinementStrategy::minStrategy() {
  return AdaptiveRefinementStrategy(
      [](std::vector<double> const &predecessorNorms, size_t numNewPoints) -> double {
        if (predecessorNorms.size() == 0)
          throw std::runtime_error("No predecessors in AdaptiveRefinementStrategy");

        double val = predecessorNorms[0];

        for (size_t i = 1; i < predecessorNorms.size(); ++i) {
          double other = predecessorNorms[i];
          if (other < val) {
            val = other;
          }
        }

        return val / static_cast<double>(numNewPoints);
      });
}

AdaptiveRefinementStrategy AdaptiveRefinementStrategy::arithmeticMeanStrategy() {
  return AdaptiveRefinementStrategy(
      [](std::vector<double> const &predecessorNorms, size_t numNewPoints) -> double {
        if (predecessorNorms.size() == 0)
          throw std::runtime_error("No predecessors in AdaptiveRefinementStrategy");

        double sum = predecessorNorms[0];

        for (size_t i = 1; i < predecessorNorms.size(); ++i) {
          sum += predecessorNorms[i];
        }

        return sum / static_cast<double>(numNewPoints);
      });
}

AdaptiveRefinementStrategy AdaptiveRefinementStrategy::geometricMeanStrategy() {
  return AdaptiveRefinementStrategy(
      [](std::vector<double> const &predecessorNorms, size_t numNewPoints) -> double {
        if (predecessorNorms.size() == 0)
          throw std::runtime_error("No predecessors in AdaptiveRefinementStrategy");

        double prod = predecessorNorms[0];

        for (size_t i = 1; i < predecessorNorms.size(); ++i) {
          prod *= predecessorNorms[i];
        }

        return pow(prod, 1.0 / static_cast<double>(numNewPoints));
      });
}

} /* namespace combigrid */
} /* namespace sgpp*/
