// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/combigrid/LevelIndexTypes.hpp>

namespace sgpp {
namespace combigrid {

/**
 * @brief a generic relevance calculator for dimensionally-adaptive combination technique
 *
 * cf. [0] Gerstner, T. and Griebel, M., 2003. Dimension–adaptive tensor–product quadrature.
 * Computing, 71(1), pp.65-87.
 * <- here, it is called "Error Estimation"
 */
class RelevanceCalculator {
 public:
  /**
   * @brief empty virtual destructor
   */
  virtual ~RelevanceCalculator() {}

  /**
   * @brief get a relevance for the subspace of LevelVector levelVector and Delta delta
   */
  virtual double calculate(const LevelVector& levelVector, double delta) const = 0;
};

}  // namespace combigrid
}  // namespace sgpp
