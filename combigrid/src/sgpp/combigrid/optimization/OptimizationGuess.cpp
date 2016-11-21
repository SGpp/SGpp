// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/optimization/OptimizationGuess.hpp>

namespace sgpp {
namespace combigrid {

bool operator<(const OptimizationGuess& first, const OptimizationGuess& second) {
  return first.fb < second.fb;
}

bool operator>(const OptimizationGuess& first, const OptimizationGuess& second) {
  return first.fb > second.fb;
}

} /* namespace combigrid */
} /* namespace sgpp */
