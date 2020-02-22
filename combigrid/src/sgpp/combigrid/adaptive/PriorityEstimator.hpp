// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/combigrid/LevelIndexTypes.hpp>

#include <map>

namespace sgpp {
namespace combigrid {

/**
 * @brief a generic priority estimator for level vectors that don't have a definite result / QoI
 * yet
 */
class PriorityEstimator {
 public:
  /**
   * @brief empty virtual destructor
   */
  virtual ~PriorityEstimator() {}

  /**
   * @brief get a priority estimate based on the downward neighbors' deltas
   *
   * @param levelVector                 the level of the level vector considered
   * @param deltasOfDownwardNeighbors   the levels and deltas of the downward neighbors of
   *                                        levelVector
   * @return double                     the priority
   */
  virtual double estimatePriority(
      const LevelVector& levelVector,
      const std::map<LevelVector, double>& deltasOfDownwardNeighbors) const = 0;
};

}  // namespace combigrid
}  // namespace sgpp
