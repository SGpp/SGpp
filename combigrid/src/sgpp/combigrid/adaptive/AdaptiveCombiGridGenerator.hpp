// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/combigrid/LevelIndexTypes.hpp>
#include <sgpp/combigrid/grid/CombinationGrid.hpp>

#include <map>
#include <vector>

namespace sgpp {
namespace combigrid {

/**
 * @brief a generic error estimator for dimensionally-adaptive combination technique
 *
 * cf. Gerstner, T. and Griebel, M., 2003. Dimension–adaptive tensor–product quadrature. Computing,
 * 71(1), pp.65-87.
 */
class ErrorEstimator {
 public:
  /**
   * @brief get an error / priority estimate for the subspace of LevelVector levelVector and Delta
   * delta
   */
  double estimate(const LevelVector& levelVector, double delta) const;
};

/**
 * @brief a generic priority calculator for subspaces that don't have a definite result yet
 */
class PriorityCalculator {
 public:
  /**
   * @brief get a priority estimate based on the downward neighbors' deltas
   *
   * @param levelVector                 the level of the subspace considered
   * @param deltasOfDownwardNeighbors   the levels and deltas of the downward neighbors of
   *                                        levelVector
   * @return double                     the priority
   */
  double calculatePriority(const LevelVector& levelVector,
                           const std::map<LevelVector, double>& deltasOfDownwardNeighbors) const;
};

/**
 * @brief The AdaptiveCombiGridGenerator is a (potentially changing) representation of a combination
 * grid that also tracks the error estimates; in contrast to a CombinationGrid, it stores the full
 * downward closed set of levels.
 *
 * cf. Gerstner, T. and Griebel, M., 2003. Dimension–adaptive tensor–product quadrature. Computing,
 * 71(1), pp.65-87.
 */
class AdaptiveCombiGridGenerator {
 public:
  /**
   * Default constructor, starting with no subspaces
   */
  AdaptiveCombiGridGenerator();

  /**
   * @brief Construct a new Adaptive Level Set object
   *
   * @param combinationGrid     start with the subspaces contained in combinationGrid
   * @param errorEstimator      an error estimator relating deltas and level vectors to an
   * error/priority estimate
   * @param priorityStrategy    a priority strategy to get the priority of a level / subspace whose
   * result we don't yet know
   */
  AdaptiveCombiGridGenerator(CombinationGrid combinationGrid, ErrorEstimator errorEstimator,
                             PriorityCalculator priorityCalculator);

  //
  /**
   * @brief Get the the currently valid combination grid consisting of the "old set"
   * (the combination grid only holds the full grid vectors with non-zero coefficients)
   */
  CombinationGrid getCombinationGrid() const;

  /**
   * @brief Get the subspacesAndResults object
   */
  const std::map<LevelVector, double>& getSubspacesAndResults() const { return subspacesAndQoI_; }

  /**
   * @brief add information / a result on a subspace of LevelVector level
   */
  void addSubspaceInfo(const LevelVector& level, double scalar);

  /**
   * @brief add the next most important subspace of known result to the old set
   *
   * @param regular   add only subspaces of the next regular level, or none
   * @return true     if there was a subspace added
   * @return false    if there was no subspace added (e.g. because all known results were in the old
   *                  set already)
   */
  bool adaptNextSubspace(bool regular = false);

  /**
   * @brief add all subspaces of known result to the old set
   */
  void adaptAllKnown();

  /**
   * @brief Get the level vectors of the active set (= admissible upward neighbors of the old set)
   *
   * @return std::vector<LevelVector> the active set
   */
  std::vector<LevelVector> getActiveSet() const;

  /**
   * @brief Get the delta belonging to a level vector
   *
   * @return double  the delta value, NaN if result is unknown
   */
  double getDelta(const LevelVector& levelVector);

  /**
   * @brief get a priority queue of elements in the active set that don't have a delta yet
   */
  std::map<LevelVector, double> getPriorityQueue() const;

  /**
   * @brief get exact value of priority / error of those elements in the active set
   * that already have a value
   */
  std::map<LevelVector, double> getErrorsOfActiveSet() const;

  /**
   * @brief get an estimate or exact value of priority of all elements in the active set
   */
  std::map<LevelVector, double> getPrioritiesOfActiveSet() const {
    auto errors = getErrorsOfActiveSet();
    auto priorityQueueAndErrors = getPriorityQueue();
    priorityQueueAndErrors.insert(errors.begin(), errors.end());
    return priorityQueueAndErrors;
  }

 private:
  // a map that holds all the levels / results
  // key: the downward-closed set of all level vectors / subspaces considered so far
  // value: the results obtained by evaluating the full grids, according to which the grid will be
  // adapted; by default initialized with NaN (in case there is no result yet)
  // TODO(pollinta) allow results to be something else than a double / scalar
  //                ( = function, vector...)
  std::map<LevelVector, double> subspacesAndQoI_;

  // the old set = the subspaces that are definitely in our combigrid already
  std::vector<LevelVector> oldSet_;

  // the error estimator used to relate delta and level vector to an "error" / importance
  ErrorEstimator errorEstimator_;  // delta, variance (level, delta -> error)

  // the priority calculator used to calculate the priority of a level whose result we don't yet
  // know (similar to errorEstimator_, but based on the deltas of the downward neighbors)
  PriorityCalculator priorityCalculator_;  // averaging, weighted
};
}  // namespace combigrid
} /* namespace sgpp */
