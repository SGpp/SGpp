// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/exception/not_implemented_exception.hpp>

#include <sgpp/combigrid/LevelIndexTypes.hpp>
#include <sgpp/combigrid/grid/CombinationGrid.hpp>

#include <algorithm>
#include <cassert>
#include <limits>
#include <list>
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
 * @brief a generic priority calculator for level vectors that don't have a definite result / QoI
 * yet
 */
class PriorityCalculator {
 public:
  /**
   * @brief get a priority estimate based on the downward neighbors' deltas
   *
   * @param levelVector                 the level of the level vector considered
   * @param deltasOfDownwardNeighbors   the levels and deltas of the downward neighbors of
   *                                        levelVector
   * @return double                     the priority
   */
  virtual double calculatePriority(
      const LevelVector& levelVector,
      const std::map<LevelVector, double>& deltasOfDownwardNeighbors) const = 0;
};

/**
 * @brief the AveragingLevelManager from @holzmudd's combigrid module
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
 * The adaptation is performed based on a scalar Quantity of Interest passed to the
 * AdaptiveCombiGridGenerator via \c addQoIInformation .
 *
 * There are -- potentially -- three categories of level vectors in an AdaptiveCombiGridGenerator,
 * those that are already selected by the adaptive algorithm ("old set"), candidates for the old set
 * which are not yet selected whose QoIs we may or may not yet know ("active set"), and others that
 * are outside either of those (whose QoI we are already storing for later use).
 *
 * When doing adaptation by \c adaptNextLevelVector , an admissible level vector is moved from the
 * active to the old set. That level vector will be selected by the \c errorEstimator based on
 * delta, a measure for how much the combined QoI will change if a particular level vector would be
 * added.
 *
 * Terminology is mostly taken from Gerstner, T. and Griebel, M., 2003. Dimension–adaptive
 * tensor–product quadrature. Computing, 71(1), pp.65-87.
 */
class AdaptiveCombiGridGenerator {
  using MapPairType = std::pair<const std::vector<unsigned int>, double>;

 public:
  /**
   * Default constructor, starting with no level vectors
   */
  AdaptiveCombiGridGenerator();

  /**
   * @brief Construct a new AdaptiveCombiGridGenerator object
   *
   * @param levelVectors           start with these level vectors, cannot be empty
   * @param errorEstimator      an error estimator relating deltas and level vectors to an
   *                                error/priority estimate
   * @param priorityStrategy    a priority strategy to get the priority of a level / subspace whose
   *                                result we don't yet know
   */
  AdaptiveCombiGridGenerator(const std::vector<LevelVector>& levelVectors,
                             std::unique_ptr<ErrorEstimator> errorEstimator,
                             std::unique_ptr<PriorityCalculator> priorityCalculator)
      : errorEstimator_(std::move(errorEstimator)),
        priorityCalculator_(std::move(priorityCalculator)) {
    assert(levelVectors.size() > 0);

    // set the minimum level vector
    auto numDimensions = levelVectors[0].size();
    minimumLevelVector_ = LevelVector(numDimensions, std::numeric_limits<level_t>::max());
    for (const auto& subspace : levelVectors) {
      for (size_t d = 0; d < numDimensions; ++d) {
        minimumLevelVector_[d] = std::min(subspace[d], minimumLevelVector_[d]);
      }
    }

    auto downwardClosedLevelSet = makeDownwardClosed(levelVectors, minimumLevelVector_);

    for (const auto& level : downwardClosedLevelSet) {
      subspacesAndQoI_[level] = std::numeric_limits<double>::quiet_NaN();
      activeSet_.push_back(level);
      adaptLevel(level);
    }
  }

  /**
   * @brief Construct a new AdaptiveCombiGridGenerator object
   *
   * @param combinationGrid     start with the subspaces contained in combinationGrid (must be at
   *                            least one)
   * @param errorEstimator      an error estimator relating deltas and level vectors to an
   *                                error/priority estimate
   * @param priorityStrategy    a priority strategy to get the priority of a level / subspace whose
   *                                result we don't yet know
   */
  static AdaptiveCombiGridGenerator fromCombinationGrid(
      const CombinationGrid& combinationGrid, std::unique_ptr<ErrorEstimator> errorEstimator,
      std::unique_ptr<PriorityCalculator> priorityCalculator) {
    std::vector<LevelVector> subspaces{};
    subspaces.resize(combinationGrid.getFullGrids().size());
    std::transform(combinationGrid.getFullGrids().begin(), combinationGrid.getFullGrids().end(),
                   subspaces.begin(),
                   [](const FullGrid& fg) -> LevelVector { return fg.getLevel(); });
    return AdaptiveCombiGridGenerator(subspaces, std::move(errorEstimator),
                                      std::move(priorityCalculator));
  }

  /**
   * @brief Get the the currently valid combination grid consisting of the "old set"
   * (the combination grid only holds the full grid vectors with non-zero coefficients)
   */
  CombinationGrid getCombinationGrid(const HeterogeneousBasis& basis) const {
    bool hasBoundary = std::find(minimumLevelVector_.begin(), minimumLevelVector_.end(), 0) !=
                       minimumLevelVector_.end();
    return CombinationGrid::fromSubspaces(oldSet_, basis, hasBoundary);
  }

  /**
   * @brief Get the subspacesAndQoIs object
   */
  const std::map<LevelVector, double>& getSubspacesAndQoIs() const { return subspacesAndQoI_; }

  /**
   * @brief add information / a result on a subspace of LevelVector level
   */
  void addQoIInformation(const LevelVector& level, double qoi) { subspacesAndQoI_[level] = qoi; }

  /**
   * @brief add the next most important subspace of known result to the old set
   *
   * @param regular   add only subspaces of the next regular level, or none
   * @return true     if there was a subspace added
   * @return false    if there was no subspace added (e.g. because all known results were in the old
   *                  set already) //TODO(pollinta): use regular parameter
   */
  bool adaptNextLevelVector(bool regular = false) {
    if (regular) {
      throw sgpp::base::not_implemented_exception("Parameter regular not yet implemented!");
    }
    auto errors = getErrorsOfActiveSet();
    auto pr = std::max_element(
        errors.begin(), errors.end(),
        [](const MapPairType& m1, const MapPairType& m2) { return m1.second < m2.second; });
    if (pr != errors.end()) {
      adaptLevel(pr->first);
      return true;
    } else {
      return false;
    }
  }

  /**
   * @brief add all subspaces of known result to the old set
   */
  bool adaptAllKnown() {
    bool levelVectorAdded = false;
    do {
      levelVectorAdded = adaptNextLevelVector();
    } while (levelVectorAdded);
    return levelVectorAdded;
  }

  /**
   * @brief Get the level vectors of the old set
   *
   * @return std::vector<LevelVector> the old set
   */
  std::vector<LevelVector> getOldSet() const { return oldSet_; }

  /**
   * @brief Get the level vectors of the active set (= admissible upward neighbors of the old set)
   *
   * @return std::list<LevelVector> the active set
   */
  std::list<LevelVector> getActiveSet() const { return activeSet_; }

  /**
   * @brief Get the minimum Level Vector object
   * 
   * @return const LevelVector& 
   */
  const LevelVector& getMinimumLevelVector() const { return minimumLevelVector_; }

  /**
   * @brief Get the delta belonging to a level vector. Delta denotes by how much the combined value
   * of QoI would be influenced if this level vector would be adapted to
   *
   * @return double  the delta value, NaN if result is unknown
   */
  double getDelta(const LevelVector& levelVector) const {

    if (subspacesAndQoI_.find(levelVector) == subspacesAndQoI_.end()) {
      return std::numeric_limits<double>::quiet_NaN();
    }

    double neighborStencilSum = 0.;
    auto levelVectorMinusOne = levelVector;
    for (size_t d = 0; d < levelVectorMinusOne.size(); ++d) {
      auto& l = levelVectorMinusOne[d];
      if (l > minimumLevelVector_[d]) {
        l -= 1;
      }
    }

    auto lowerHypercube = hyperCubeOfLevelVectors(levelVector, levelVectorMinusOne);
    lowerHypercube.pop_back();
    // TODO(pollinta): simplify Hamming distance calculation
    for (size_t i = 0; i < lowerHypercube.size(); ++i) {
      auto hypercubeElement = subspacesAndQoI_.find(lowerHypercube[i]);
      assert(hypercubeElement != subspacesAndQoI_.end());
      auto hammingDistance = 0;
      for (size_t d = 0; d < levelVector.size(); ++d) {
        hammingDistance += levelVector[d] - lowerHypercube[i][d];
      }
      neighborStencilSum += subspacesAndQoI_.at(lowerHypercube[i]) * std::pow(-1, hammingDistance);
    }
    return subspacesAndQoI_.at(levelVector) - neighborStencilSum;
  }

  /**
   * @brief get a priority queue of elements in the active set that don't have a result / QoI /
   * delta yet
   */
  std::map<LevelVector, double> getPriorityQueue() const;

  /**
   * @brief get exact value of priority / error of those elements in the active set
   * that already have a QoI value
   */
  std::map<LevelVector, double> getErrorsOfActiveSet() const {
    std::map<LevelVector, double> errors{};
    for (const auto& levelVector : activeSet_) {
      const auto delta = getDelta(levelVector);
      if (!std::isnan(delta)) {
        errors[levelVector] = errorEstimator_->estimate(levelVector, delta);
      }
    }
    return errors;
  }

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
  /**
   * @brief a level vector is admissible, if all potential lower neighbors are already in the old
   * set
   */
  bool isAdmissible(const LevelVector& level) const {
    for (size_t d = 0; d < minimumLevelVector_.size(); ++d) {
      if (level[d] > minimumLevelVector_[d]) {
        auto neighborLevel = level;
        neighborLevel[d] -= 1;
        bool isInOldSet = std::find(oldSet_.begin(), oldSet_.end(), neighborLevel) != oldSet_.end();
        if (!isInOldSet) {
          return false;
        }
      }
    }
    return true;
  }

  /**
   * @brief add forward neighbors of \c level to active set, if admissible
   */
  void addNeighborsToActiveSet(const LevelVector& level) {
    for (size_t d = 0; d < level.size(); ++d) {
      auto neighborLevel = level;
      neighborLevel[d] += 1;
      if (isAdmissible(neighborLevel)) {
        activeSet_.push_back(neighborLevel);
      }
    }
  }

  /**
   * @brief adapt to level vector \c level = move from active to old set
   */
  void adaptLevel(const LevelVector& level) {
    assert(std::find(oldSet_.begin(), oldSet_.end(), level) == oldSet_.end());
    assert(std::find(activeSet_.begin(), activeSet_.end(), level) != activeSet_.end());

    oldSet_.push_back(level);
    activeSet_.remove(level);
    addNeighborsToActiveSet(level);
  }

  // a map that holds all the levels / results
  // key: the downward-closed set of all level vectors / subspaces considered so far
  // value: the results obtained by evaluating the full grids, according to which the grid will be
  // adapted; by default initialized with NaN (in case there is no result yet)
  // TODO(pollinta) allow results to be something else than a double / scalar ... how to do this?
  //                ( = function, vector...)
  std::map<LevelVector, double> subspacesAndQoI_;

  // the minimum level vector, from where the adaptation is started
  LevelVector minimumLevelVector_;

  // the old set = the level vectors that are definitely in our combigrid already
  std::vector<LevelVector> oldSet_;

  // the active set = the level vectors that may be added to our combigrid next
  std::list<LevelVector> activeSet_;

  // the error estimator used to relate delta and level vector to an "error" / importance
  std::unique_ptr<ErrorEstimator> errorEstimator_;  // delta, variance (level, delta -> error)

  // the priority calculator used to calculate the priority of a level whose result we don't yet
  // know (similar to errorEstimator_, but based on the deltas of the downward neighbors)
  std::unique_ptr<PriorityCalculator> priorityCalculator_;  // averaging, weighted
};
}  // namespace combigrid
} /* namespace sgpp */
