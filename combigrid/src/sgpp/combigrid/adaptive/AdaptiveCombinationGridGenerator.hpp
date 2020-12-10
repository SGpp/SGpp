// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/combigrid/LevelIndexTypes.hpp>
#include <sgpp/combigrid/adaptive/AveragingPriorityEstimator.hpp>
#include <sgpp/combigrid/adaptive/DataAbstraction.hpp>
#include <sgpp/combigrid/adaptive/PriorityEstimator.hpp>
#include <sgpp/combigrid/adaptive/RelevanceCalculator.hpp>
#include <sgpp/combigrid/adaptive/WeightedRelevanceCalculator.hpp>
#include <sgpp/combigrid/grid/CombinationGrid.hpp>

#include <functional>
#include <limits>
#include <list>
#include <map>
#include <memory>
#include <utility>
#include <vector>

// template<>
// class T std::plus<class T>(const & T a, const & T b){
//     T result;
//     for (i=0; i< a.size(); ++i){
//         result[i]= a[i]+b[i];
//     }
//     return result;
// }

namespace sgpp {
namespace combigrid {

/**
 * @brief The AdaptiveCombinationGridGenerator is a (potentially changing) representation of a
 * combination grid that also tracks the Quantities of Interest; in contrast to a CombinationGrid,
 * it stores the full downward closed set of levels.
 *
 * The adaptation is performed based on a scalar Quantity of Interest passed to the
 * AdaptiveCombinationGridGenerator via \c setQoIInformation .
 *
 * There are -- potentially -- three categories of level vectors in an
 * AdaptiveCombinationGridGenerator, those that are already selected by the adaptive algorithm ("old
 * set"), candidates for the old set which are not yet selected whose QoIs we may or may not yet
 * know ("active set"), and others that are outside either of those (whose QoI we are already
 * storing for later use).
 *
 * When doing adaptation by \c adaptNextLevelVector , an admissible level vector is moved from the
 * active to the old set. That level vector will be selected by the \c relevanceCalculator based on
 * delta, a measure for how much the combined QoI will change if a particular level vector would be
 * added.
 *
 * If grid evaluations are expensive and one would like to know which QoIs should be computed next,
 * a priority queue can be obtained by \c getPriorityQueue , which uses the \c priorityEstimator to
 * infer a priority for the active set levels from the QoIs of the downward neighbors.
 *
 * Terminology is mostly taken from Gerstner, T. and Griebel, M., 2003. Dimension–adaptive
 * tensor–product quadrature. Computing, 71(1), pp.65-87.
 */
template <typename T>
class AdaptiveGenerator {
  typedef std::pair<const std::vector<unsigned int>, double> MapPairType;

 public:
  /**
   * @brief Construct a new AdaptiveCombinationGridGenerator object
   *
   * @param levelVectors           start with these level vectors, cannot be empty
   * @param QoIValues              the QoI values corresponding to levelVectors, same length as
   *                                 levelVectors
   * @param summationFunction      the summation function by which results are combined
   * @param relevanceCalculator a relevance calculator relating deltas and level vectors to an
   *                                "error"/relevance estimate
   * @param priorityEstimator   a priority estimator to get the priority of a level / subspace whose
   *                                result we don't yet know (similar to relevanceCalculator, but
   *                                based on the deltas of the downward neighbors instead of the
   *                                level's own delta)
   */
  AdaptiveGenerator(const std::vector<LevelVector>& levelVectors, const std::vector<T>&& QoIValues,
                    std::function<T(T, T)> summationFunction,
                    std::shared_ptr<RelevanceCalculator<T>> relevanceCalculator =
                        std::shared_ptr<RelevanceCalculator<T>>(new WeightedRelevanceCalculator()),
                    std::shared_ptr<PriorityEstimator<T>> priorityEstimator =
                        std::shared_ptr<PriorityEstimator<T>>(new AveragingPriorityEstimator()));

  /**
   * like above, with default QoI values if they are not supplied
   */
  AdaptiveGenerator(const std::vector<LevelVector>& levelVectors,
                    std::function<T(T, T)> summationFunction,
                    std::shared_ptr<RelevanceCalculator<T>> relevanceCalculator =
                        std::shared_ptr<RelevanceCalculator<T>>(new WeightedRelevanceCalculator()),
                    std::shared_ptr<PriorityEstimator<T>> priorityEstimator =
                        std::shared_ptr<PriorityEstimator<T>>(new AveragingPriorityEstimator()))
      : AdaptiveGenerator(levelVectors, std::vector<T>(levelVectors.size(), getQuietNan<T>()),
                          summationFunction, relevanceCalculator, priorityEstimator) {}

  /**
   * like above, but setting the summationFunction to std::plus<T>() by default
   */
  AdaptiveGenerator(const std::vector<LevelVector>& levelVectors,
                    std::shared_ptr<RelevanceCalculator<T>> relevanceCalculator =
                        std::shared_ptr<RelevanceCalculator<T>>(new WeightedRelevanceCalculator()),
                    std::shared_ptr<PriorityEstimator<T>> priorityEstimator =
                        std::shared_ptr<PriorityEstimator<T>>(new AveragingPriorityEstimator()))
      : AdaptiveGenerator(levelVectors, std::plus<T>(), relevanceCalculator, priorityEstimator) {}

  /**
   * @brief Construct a new AdaptiveGenerator object
   *
   * @param combinationGrid     start with the subspaces contained in combinationGrid (must be at
   *                            least one)
   * @param QoIValues           the QoI values corresponding to levelVectors, same length as
   *                              the number of FullGrids in combinationGrid
   * @param summationFunction   the summation function by which results are combined
   * @param relevanceCalculator a relevance calculator relating deltas and level vectors to an
   *                                "error"/relevance estimate
   * @param priorityEstimator   a priority estimator to get the priority of a level / subspace whose
   *                                result we don't yet know (similar to relevanceCalculator, but
   *                                based on the deltas of the downward neighbors instead of the
   *                                level's own delta)
   */
  static AdaptiveGenerator fromCombinationGrid(
      const CombinationGrid& combinationGrid, const std::vector<T>&& QoIValues,
      std::function<T(T, T)> summationFunction = std::plus<T>(),
      std::shared_ptr<RelevanceCalculator<T>> relevanceCalculator =
          std::shared_ptr<RelevanceCalculator<T>>(new WeightedRelevanceCalculator()),
      std::shared_ptr<PriorityEstimator<T>> priorityEstimator =
          std::shared_ptr<PriorityEstimator<T>>(new AveragingPriorityEstimator()));

  /**
   * like above, with default QoI values if they are not supplied
   */
  static AdaptiveGenerator fromCombinationGrid(
      const CombinationGrid& combinationGrid,
      std::function<T(T, T)> summationFunction = std::plus<T>(),
      std::shared_ptr<RelevanceCalculator<T>> relevanceCalculator =
          std::shared_ptr<RelevanceCalculator<T>>(new WeightedRelevanceCalculator()),
      std::shared_ptr<PriorityEstimator<T>> priorityEstimator =
          std::shared_ptr<PriorityEstimator<T>>(new AveragingPriorityEstimator())) {
    return fromCombinationGrid(
        combinationGrid, std::vector<T>(combinationGrid.getFullGrids().size(), getQuietNan<T>()),
        summationFunction, relevanceCalculator, priorityEstimator);
  }

  /**
   * @brief Get the the currently valid combination grid consisting of the "old set"
   * (the combination grid only holds the full grid vectors with non-zero coefficients)
   */
  CombinationGrid getCombinationGrid(const HeterogeneousBasis& basis,
                                     bool hasBoundary = true) const;

  /**
   * @brief Get the subspacesAndQoIs object
   */
  const std::map<LevelVector, T>& getSubspacesAndQoIs() const { return subspacesAndQoI; }

  /**
   * @brief set QoI information / a result for LevelVector level
   */
  void setQoIInformation(const LevelVector& level, T qoi) { subspacesAndQoI[level] = qoi; }

  /**
   * @brief add the next most important subspace of known result to the old set
   *
   * @param regular   add only subspaces of the next regular level, or none
   * @return true     if there was a subspace added
   * @return false    if there was no subspace added (e.g. because all known results were in the old
   *                  set already) //TODO(pollinta): use regular parameter
   */
  bool adaptNextLevelVector(bool regular = false);

  /**
   * @brief add all subspaces of known result to the old set
   */
  bool adaptAllKnown();

  /**
   * @brief Get the level vectors of the old set
   *
   * @return std::vector<LevelVector> the old set
   */
  std::vector<LevelVector> getOldSet() const { return oldSet; }

  /**
   * @brief Get current result based on the old set
   *
   * @return the combined QoI
   */
  T getCurrentResult() const;

  /**
   * @brief Get the level vectors of the active set (= admissible upward neighbors of the old set)
   *
   * @return std::vector<LevelVector> the active set
   */
  std::vector<LevelVector> getActiveSet() const;

  /**
   * @brief Get the level vectors of the active set (= admissible upward neighbors of the old set)
   *
   * @return std::vector<LevelVector> all levels (old and active)
   */
  std::vector<LevelVector> getLevels() const;

  /**
   * @brief Get the minimum Level Vector object
   *
   * @return const LevelVector&
   */
  const LevelVector& getMinimumLevelVector() const { return minimumLevelVector; }

  /**
   * @brief Get the delta belonging to a level vector. Delta denotes by how much the combined value
   * of QoI would be influenced if this level vector would be adapted to
   *
   * @return T  the delta value, NaN if result is unknown
   */
  T getDelta(const LevelVector& levelVector) const;

  /**
   * @brief Get the deltas belonging to the list of supplied level vectors
   *
   * @return T  the delta value, NaN if result is unknown
   */
  std::vector<T> getDeltas(const std::vector<LevelVector>& levelVectors) const;

  /**
   * @brief get the levels and priority of elements in the active set that don't have a result / QoI
   * value / delta yet
   */
  std::map<LevelVector, double> getPriorities() const;

  /**
   * @brief get a priority queue of elements in the active set that don't have a result / QoI /
   * delta yet
   */
  std::vector<LevelVector> getPriorityQueue() const;

  /**
   * @brief get exact value of relevance / "error" of those elements in the active set
   * that already have a QoI value
   */
  std::map<LevelVector, double> getRelevanceOfActiveSet() const;

  /**
   * @brief adapt to level vector \c level = move from active to old set
   */
  void adaptLevel(const LevelVector& level);

 private:
  /**
   * @brief a level vector is admissible, if all potential lower neighbors are already in the old
   * set
   */
  bool isAdmissible(const LevelVector& level) const;

  /**
   * @brief add forward neighbors of \c level to active set, if admissible
   */
  void addNeighborsToActiveSet(const LevelVector& level);

  // a map that holds all the levels / results
  // key: the downward-closed set of all level vectors / subspaces considered so far
  // value: the results obtained by evaluating the full grids, according to which the grid will be
  // adapted; by default initialized with NaN (in case there is no result yet)
  // TODO(pollinta) allow results to be something else than a double / scalar ... how to do this?
  //                ( = function, vector...)
  std::map<LevelVector, T> subspacesAndQoI;

  // the function by which the individual QoIs will be combined to compute deltas and to give the
  // total result. Default is linear summation, variance adaptation would use quadratic summation
  std::function<T(T, T)> summationFunction;

  // the minimum level vector, from where the adaptation is started
  LevelVector minimumLevelVector;

  // the old set = the level vectors that are definitely in our combigrid already
  std::vector<LevelVector> oldSet;

  // the active set = the level vectors that may be added to our combigrid next
  std::list<LevelVector> activeSet;

  // the relevance calculator used to relate delta and level vector to an "error" / relevance
  std::shared_ptr<RelevanceCalculator<T>> relevanceCalculator;

  // the priority estimator used to estimate the priority of a level whose result we don't yet
  // know (similar to relevanceCalculator, but based on the deltas of the downward neighbors
  // instead of the level's own delta)
  std::shared_ptr<PriorityEstimator<T>> priorityEstimator;  // averaging, weighted
};

typedef AdaptiveGenerator<double> ScalarAdaptiveGenerator;
typedef ScalarAdaptiveGenerator AdaptiveCombinationGridGenerator;
}  // namespace combigrid
}  // namespace sgpp
