// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/combigrid/LevelIndexTypes.hpp>
#include <sgpp/combigrid/grid/CombinationGrid.hpp>

#include <algorithm>
#include <cassert>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <utility>
#include <vector>

namespace sgpp {
namespace combigrid {

/**
 * @brief a generic error estimator for dimensionally-adaptive combination technique
 *
 * cf. [0] Gerstner, T. and Griebel, M., 2003. Dimension–adaptive tensor–product quadrature.
 * Computing, 71(1), pp.65-87.
 */
class ErrorEstimator {
 public:
  /**
   * @brief get an error / priority estimate for the subspace of LevelVector levelVector and Delta
   * delta
   */
  virtual double estimate(const LevelVector& levelVector, double delta) const = 0;
};

/**
 * @brief the weighted error estimator introduced in [0]
 */
class WeightedErrorEstimator : public ErrorEstimator {
 public:
  explicit WeightedErrorEstimator(double weightErrorInRelationToNumberOfPoints = 1.)
      : weightErrorInRelationToNumberOfPoints_(weightErrorInRelationToNumberOfPoints) {}

  double estimate(const LevelVector& levelVector, double delta) const override {
    auto numPoints = static_cast<index_t>(1)
                     << std::accumulate(levelVector.begin(), levelVector.end(), 0);
    return std::max(weightErrorInRelationToNumberOfPoints_ * delta,
                    (1 - weightErrorInRelationToNumberOfPoints_) / static_cast<double>(numPoints));
  }

 private:
  double weightErrorInRelationToNumberOfPoints_;
};

/**
 * @brief a generic priority calculator for subspaces that don't have a definite result / QoI yet
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
  virtual double calculatePriority(
      const LevelVector& levelVector,
      const std::map<LevelVector, double>& deltasOfDownwardNeighbors) const = 0;
};

/**
 * @brief the AveragingLevelManager from holzmudd's combigrid module
 */
class AveragingPriorityCalculator : public PriorityCalculator {
 public:
  double calculatePriority(const LevelVector& levelVector,
                           const std::map<LevelVector, double>& deltasOfDownwardNeighbors) const {
    auto normDividedByNumberOfPoints =
        [](double& accumulateResult,
           const std::pair<const std::vector<unsigned int>, double>& mapEntry) {
          auto sumOfLevelVector = static_cast<index_t>(1) << std::accumulate(
                                      mapEntry.first.begin(), mapEntry.first.end(), 0);
          accumulateResult += mapEntry.second / sumOfLevelVector;
          return accumulateResult;
        };
    auto sumOfNormDividedByNumberOfPoints =
        std::accumulate(deltasOfDownwardNeighbors.begin(), deltasOfDownwardNeighbors.end(), 0.,
                        normDividedByNumberOfPoints);
    return sumOfNormDividedByNumberOfPoints / static_cast<double>(deltasOfDownwardNeighbors.size());
  }
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
   * @param subspaces           start with these subspaces, cannot be empty
   * @param errorEstimator      an error estimator relating deltas and level vectors to an
   *                                error/priority estimate
   * @param priorityStrategy    a priority strategy to get the priority of a level / subspace whose
   *                                result we don't yet know
   * @param QoIs                the results / QoIs, if already known for all subspaces
   */
  AdaptiveCombiGridGenerator(const std::vector<LevelVector>& subspaces,
                             std::unique_ptr<ErrorEstimator> errorEstimator,
                             std::unique_ptr<PriorityCalculator> priorityCalculator,
                             std::vector<double> QoIs = std::vector<double>())
      : errorEstimator_(std::move(errorEstimator)),
        priorityCalculator_(std::move(priorityCalculator)) {
    assert((QoIs.size() == 0) || (QoIs.size() == subspaces.size()));
    assert(subspaces.size() > 0);

    if (QoIs.size() == 0) {
      for (const auto& subspace : subspaces) {
        subspacesAndQoI_[subspace] = std::numeric_limits<double>::quiet_NaN();
      }
    } else {
      for (size_t i = 0; i < subspaces.size(); ++i) {
        subspacesAndQoI_[subspaces[i]] = QoIs[i];
      }
    }
    // set the minimum level vector
    auto numDimensions = subspaces[0].size();
    minimumLevelVector_ = LevelVector(numDimensions, std::numeric_limits<level_t>::max());
    for (const auto& subspace : subspaces) {
      for (size_t d = 0; d < numDimensions; ++d) {
        minimumLevelVector_[d] = std::min(subspace[d], minimumLevelVector_[d]);
      }
    }
  }

  /**
   * @brief Construct a new Adaptive Level Set object
   *
   * @param combinationGrid     start with the subspaces contained in combinationGrid (must be at
   *                            least one)
   * @param errorEstimator      an error estimator relating deltas and level vectors to an
   *                                error/priority estimate
   * @param priorityStrategy    a priority strategy to get the priority of a level / subspace whose
   *                                result we don't yet know
   * @param QoIs                the results / QoIs, if already known for all subspaces in
   * combinationGrid
   */
  static AdaptiveCombiGridGenerator fromCombinationGrid(
      const CombinationGrid& combinationGrid, std::unique_ptr<ErrorEstimator> errorEstimator,
      std::unique_ptr<PriorityCalculator> priorityCalculator,
      std::vector<double> QoIs = std::vector<double>()) {
    std::vector<LevelVector> subspaces{};
    std::transform(combinationGrid.getFullGrids().begin(), combinationGrid.getFullGrids().end(),
                   subspaces.begin(),
                   [](const FullGrid& fg) -> LevelVector { return fg.getLevel(); });
    return AdaptiveCombiGridGenerator(subspaces, std::move(errorEstimator),
                                      std::move(priorityCalculator), QoIs);
  }

  /**
   * @brief Get the the currently valid combination grid consisting of the "old set"
   * (the combination grid only holds the full grid vectors with non-zero coefficients)
   */
  CombinationGrid getCombinationGrid() const;

  /**
   * @brief Get the subspacesAndQoIs object
   */
  const std::map<LevelVector, double>& getSubspacesAndQoIs() const { return subspacesAndQoI_; }

  /**
   * @brief add information / a result on a subspace of LevelVector level
   */
  void addSubspaceInfo(const LevelVector& level, double qoi) { subspacesAndQoI_[level] = qoi; }

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
  double getDelta(const LevelVector& levelVector) {
    using sgpp::base::operator<<;

    double neighborStencilSum = 0.;
    auto levelVectorMinusOne = levelVector;
    for (auto& l : levelVectorMinusOne) {
      l -= 1;  // TODO < 0 or minimum?
    }

    auto lowerHypercube = hyperCubeOfLevelVectors(levelVector, levelVectorMinusOne);
    // todo remove levelvector from here
    for (size_t i = 0; i < lowerHypercube.size(); ++i) {
      auto hypercubeElement = subspacesAndQoI_.find(lowerHypercube[i]);
      assert(hypercubeElement != subspacesAndQoI_.end());
      auto hammingDistance = 0;
      for (size_t d = 0; d < levelVector.size(); ++d) {
        hammingDistance += levelVector[d] - lowerHypercube[i][d];
      }
      neighborStencilSum += subspacesAndQoI_[lowerHypercube[i]] * std::pow(-1, hammingDistance);
    }
    return subspacesAndQoI_[levelVector] - neighborStencilSum;
  }

  /**
   * @brief get a priority queue of elements in the active set that don't have a result / QoI /
   * delta yet
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
  // TODO(pollinta) allow results to be something else than a double / scalar ... how to do this?
  //                ( = function, vector...)
  std::map<LevelVector, double> subspacesAndQoI_;

  // the minimum level vector, from where the adaptation is started
  LevelVector minimumLevelVector_;

  // the old set = the subspaces that are definitely in our combigrid already
  std::vector<LevelVector> oldSet_;

  // the error estimator used to relate delta and level vector to an "error" / importance
  std::unique_ptr<ErrorEstimator> errorEstimator_;  // delta, variance (level, delta -> error)

  // the priority calculator used to calculate the priority of a level whose result we don't yet
  // know (similar to errorEstimator_, but based on the deltas of the downward neighbors)
  std::unique_ptr<PriorityCalculator> priorityCalculator_;  // averaging, weighted
};
}  // namespace combigrid
} /* namespace sgpp */
