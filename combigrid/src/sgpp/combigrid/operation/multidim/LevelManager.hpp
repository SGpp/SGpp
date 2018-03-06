// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>

#include <sgpp/combigrid/common/BoundedSumMultiIndexIterator.hpp>
#include <sgpp/combigrid/definitions.hpp>
#include <sgpp/combigrid/operation/multidim/AbstractLevelEvaluator.hpp>
#include <sgpp/combigrid/operation/multidim/AdaptiveRefinementStrategy.hpp>
#include <sgpp/combigrid/operation/multidim/LevelHelpers.hpp>
#include <sgpp/combigrid/serialization/TreeStorageSerializationStrategy.hpp>
#include <sgpp/combigrid/storage/AbstractMultiStorage.hpp>
#include <sgpp/combigrid/storage/tree/TreeStorage.hpp>

#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <queue>
#include <string>
#include <unordered_set>
#include <vector>

namespace sgpp {
namespace combigrid {

/**
 * This class provides functionality to conveniently add levels to a CombigridEvaluator.
 * It implements most of the functionality, however, the adaptive level generation strategy has to
 * be specified by deriving from this class and overriding the protected method computePriority().
 *
 * Via any LevelManager subclass, one can add levels adaptively or with a regular structure (bounded
 * 1-norm) to a CombigridEvaluator, and one can even choose to use multiple threads to precompute
 * the function values. Note that the parallelization only concerns the function evaluations and not
 * the other calculations. Thus, if your function can be computed in a very short amount of time
 * (for example (x, y) |-> exp(x) * sin(y)), the parallel evaluation might be slower than the normal
 * evaluation. Of course, if you choose to use parallel evaluation, you have to ensure that your
 * function supports multiple calls in parallel.
 *
 * The parallel evaluation functionality is designed to use all available computing power. In
 * particular, parallel adaptive level generation may also start computing levels before their
 * predecessors are finished, even though value of the level cannot be estimated as precisely as if
 * the predecessors were finished.
 */
class LevelManager {
 protected:
  // data structures for adaptive refinement
  /**
   * Priority queue ranking levels that are not yet evaluated but ready for evaluation by their
   * priority.
   */
  MultiIndexQueue queue;

  /**
   * this object stores information on the levels during adaptive refinement
   */
  std::shared_ptr<LevelInfos> infoOnAddedLevels;

  /**
   * Stores level information for already visited (!= computed) levels.
   */
  std::shared_ptr<TreeStorage<std::shared_ptr<LevelInfo>>> levelData;

  /**
   * Dimensionality of the problem.
   */
  size_t numDimensions;

  /**
   * CombigridEvaluator (hidden behind an interface to avoid carrying around its template parameter
   * here).
   */
  std::shared_ptr<AbstractLevelEvaluator> combiEval;

  /**
   * Mutex that is shared with all involved objects for evaluation parts that require mutual
   * exclusion.
   */
  std::shared_ptr<std::recursive_mutex> managerMutex;

  /**
   * Defines if statistics on the refinement process are collected or not.
   */
  bool collectStats;

  /**
   * By implementing this method in a derived class, the adaption can be customized.
   * It should yield a priority for the given level. Levels with higher priority value are added
   * earlier.
   */
  virtual double computePriority(MultiIndex const &level) = 0;

  /**
   * Initializes the data structures for adaptive level generation
   */
  virtual void initAdaption();

  /**
   * Adds successors of the given level to levelData if all of their other predecessors are also in
   * levelData.
   */
  virtual void tryAddSuccessors(MultiIndex const &level);

  /**
   * Adds a level to levelData if all of its predecessors are already in levelData.
   */
  virtual void tryAddLevel(MultiIndex const &level);

  /**
   * Puts the given level into the priority queue (uses computePriority() first to compute its
   * priority).
   */
  virtual void addToQueue(MultiIndex const &level, std::shared_ptr<LevelInfo> levelInfo);

  /**
   * Helper function returning all predecessor level multi-indices. Example: For level = (1, 0, 2),
   * it would return (0, 0, 2) and (1, 0, 1).
   */
  virtual std::vector<MultiIndex> getPredecessors(MultiIndex const &level);

  /**
   * Helper function returning all successor level multi-indices. Example: For level = (1, 0, 2), it
   * would return (2, 0, 2), (1, 1, 2) and (1, 0, 3).
   */
  virtual std::vector<MultiIndex> getSuccessors(MultiIndex const &level);

  /**
   * This function should be called between getting a level from the top of the priority queue and
   * starting computations on it. This function adapts the status of the level and its successors
   * accordingly.
   */
  virtual void beforeComputation(MultiIndex const &level);

  /**
   * This function should be called after all desired precomputations (i.e. none in the
   * single-threaded setting) for a level are done.
   * It adjusts the statuses of this level and its successors. Furthermore, it triggers calling
   * addLevel() (via predecessorsCompleted()) if all predecessor levels have already been added. If
   * this is not the case, then the level will be added via predecessorsCompleted() when the last
   * missing predecessor completes.
   */
  virtual void afterComputation(MultiIndex const &level);

  /**
   * This function is called after all desired precomputations (i.e. none in the single-threaded
   * setting) are done AND all the predecessor levels have been added. It calls addLevel() for the
   * given level, adjusts the status in levelData and calls predecessorsCompleted() on the
   * successors if their precomputations are already done. If the precomputations are not done yet,
   * it adjusts their priority using updatePriority() because there is now more information
   * available.
   */
  virtual void predecessorsCompleted(MultiIndex const &level);

  /**
   * Updates the priority of a given level in the priority queue based on the currently available
   * information.
   */
  virtual void updatePriority(MultiIndex const &level, std::shared_ptr<LevelInfo> levelInfo);

  /**
   * @return a set of level multi-indices. The levels are enumerated with increasing 1-norm until
   * the total number of necessary function evaluations would exceed the given limit maxNumPoints.
   * For example, this might return the levels (0, 0), (1, 0), (0, 1), (2, 0), (1, 1) and omit the
   * level (0, 2) because that would mean more function evaluations than specified.
   *
   * All functions with a point bound here generate levels until that point bound would be exceeded.
   * None of the functions attempts to add levels with few points to reach the point bound as close
   * as possible!
   */
  std::vector<MultiIndex> getRegularLevelsByNumPoints(size_t maxNumPoints);

  /**
   * Does precomputations for all given levels in parallel with the specified number of threads.
   * This can only be used if the level structure is known beforehand, e. g. for regular level
   * structures.
   */
  void precomputeLevelsParallel(std::vector<MultiIndex> const &levels, size_t numThreads);

  /**
   * Adds all the given levels.
   */
  void addLevels(std::vector<MultiIndex> const &levels);

  /**
   * Adds the level to the combigrid evaluator
   * @param level level
   */
  bool addLevelToCombiEval(const MultiIndex &level);

  /**
   * update statistics on added levels
   */
  void updateStats();

  /**
   * update stats for levels that have been added in the current iteration
   * @param level level to be added to the stats
   */
  void addStats(const MultiIndex &level);

 public:
  /**
   * Constructor. The CombigridEvaluator (or another derived class of AbstractLevelEvaluator) has to
   * be passed.
   * @param levelEvaluator combigrid evaluator
   * @param collectStats sets if the collection of statistics should be enabled or not during
   * refinement
   */
  explicit LevelManager(std::shared_ptr<AbstractLevelEvaluator> levelEvaluator,
                        bool collectStats = true);

  /**
   * Default constructor. If this is used, setLevelEvaluator() has to be called before adding any
   * levels.
   */
  LevelManager();

  virtual ~LevelManager();

  virtual std::shared_ptr<LevelManager> clone() = 0;

  /**
   * Sets the level evaluator (normally a CombigridEvaluator).
   */
  void setLevelEvaluator(std::shared_ptr<AbstractLevelEvaluator> levelEvaluator);

  /**
   * @param q: Maximum 1-norm of the level-multi-index, where the levels start from 0 (not from 1 as
   * in most papers).
   * If you have a norm w with levels starting from 1, simply use q = w - numDim().
   */
  void addRegularLevels(size_t q);

  /**
   * Adds a set of level multi-indices. The levels are enumerated with increasing 1-norm until
   * the total number of necessary function evaluations would exceed the given limit maxNumPoints.
   * Levels which are already computed are not counted in the number of points, so if many levels
   * have already been computed, the counting starts from their successors.
   * For example, this might add the levels (0, 0), (1, 0), (0, 1), (2, 0), (1, 1) and omit the
   * level (0, 2) because that would mean more function evaluations than specified.
   *
   * All functions with a point bound here generate levels until that point bound would be exceeded.
   * None of the functions attempts to add levels with few points to reach the point bound as close
   * as possible!
   */
  void addRegularLevelsByNumPoints(size_t maxNumPoints);

  /**
   * Does the same as addRegularLevels(), but with parallel precomputation of function values.
   * @param q  Maximum 1-norm of the level-multi-index, where the levels start from 0 (not from 1 as
   * in most papers).
   * @param numThreads number of threads that should be used for computation
   * If you have a norm w with levels starting from 1, simply use q = w - numDims().
   */
  void addRegularLevelsParallel(size_t q, size_t numThreads);

  /**
   * Does the same as addRegularLevelsByNumPoints(), but with parallel precomputation of function
   * values.
   */
  void addRegularLevelsByNumPointsParallel(size_t maxNumPoints, size_t numThreads);

  /**
   * @return the dimensionality of the problem.
   */
  virtual size_t numDims() const;

  /**
   * @return a TreeStorage which contains an entry at an index i iff the level i has been added to
   * this CombigridEvaluator.
   */
  std::shared_ptr<TreeStorage<uint8_t>> getLevelStructure() const;

  /**
   * @param q: Maximum 1-norm of the level-multi-index, where the levels start from 0 (not from 1 as
   * in most papers).
   * If you have a norm w with levels starting from 1, simply use q = w - numDims().
   */
  std::vector<MultiIndex> getRegularLevels(size_t q);

  /**
   * @return the serialized version of getLevelStructure()
   */
  std::string getSerializedLevelStructure() const;

  /**
   * Adds all levels for which an entry is contained in storage. "Inverse" operation to
   * getLevelStructure().
   */
  void addLevelsFromStructure(std::shared_ptr<TreeStorage<uint8_t>> storage);

  /**
   * Does the same as addLevelsFromStructure(), but with parallel precomputation of function values
   * using numThreads threads.
   */
  void addLevelsFromStructureParallel(std::shared_ptr<TreeStorage<uint8_t>> storage,
                                      size_t numThreads = 4);

  /**
   * Equivalent to deserializing serializedStructure and then calling addLevelsFromStructure().
   * "Inverse" operation to getSerializedLevelStructure().
   */
  void addLevelsFromSerializedStructure(std::string serializedStructure);

  /**
   * Does the same as addLevelsFromSerializedStructure(), but with parallel precomputation of
   * function values using numThreads threads.
   */
  void addLevelsFromSerializedStructureParallel(std::string serializedStructure, size_t numThreads);

  /**
   * Adds levels in an adaptive manner, such that the given maximum number of function evaluations
   * (grid points) is not exceeded. The adaption strategy depends on the particular implementation
   * of the priority function provided by the subclass of LevelManager.
   *
   * All functions with a point bound here generate levels until that point bound would be exceeded.
   * None of the functions attempts to add levels with few points to reach the point bound as close
   * as possible!
   */
  virtual void addLevelsAdaptive(size_t maxNumPoints);

  /**
   * Does the same as addLevelsAdaptive(), but with parallel function evaluations.
   */
  virtual void addLevelsAdaptiveParallel(size_t maxNumPoints, size_t numThreads);

  /**
   * @return a vector with all grid points where the function has been evaluated (without
   * duplicates).
   */
  std::vector<base::DataVector> getAllGridPoints() { return combiEval->getAllGridPoints(); }

  /**
   * @return a DataMatrix with all grid points where the function has been evaluated in its columns
   * (without duplicates)
   */
  base::DataMatrix getGridPointMatrix() { return combiEval->getGridPointMatrix(); }

  /**
   * @return the total number of different (multi-dimensional) grid points that have been used for
   * the current evaluation. This number is reset when clear() is called on the CombigridEvaluator.
   * This method is currently not optimized and can be slow!
   */
  size_t numGridPoints() { return combiEval->getAllGridPoints().size(); }

  /**
   * @return the maximum number of grid points that have to be evaluated for a regular grid with
   * maximum multi-index L1 norm of q. For nested grids, this estimate should be exact.
   */
  size_t maxNumPointsForRegular(size_t q) { return combiEval->maxNumPointsForRegular(q); }

  /**
   * Calls addLevel() on the underlying CombigridEvaluator.
   */
  void addLevel(MultiIndex const &level) { addLevelToCombiEval(level); }

  /**
   * Queue based addLevel-type function
   */
  void addLevelsAdaptiveByNumLevels(size_t numLevels = 1);

  /**
   * writes a given level structure to a matrix
   * @param levelstructure the level structure
   * @param numDims number of dimensions
   * @return matrix containing the level stucture
   */
  sgpp::base::DataMatrix convertLevelStructureToMatrix(
      std::shared_ptr<sgpp::combigrid::TreeStorage<uint8_t>> const &levelstructure, size_t numDims);

  /**
   * prints a given level structure
   * @param levelstructure the level stucture
   */
  void printLevelStructure(
      std::shared_ptr<sgpp::combigrid::TreeStorage<uint8_t>> const &levelstructure);

  /**
   * @return An upper bound for the number of points (function evaluations) used for the current
   * computation. This bound is exact if nesting is used or if otherwise each grid point only occurs
   * in exactly one level.
   */
  size_t getUpperPointBound() const { return combiEval->getUpperPointBound(); }

  /**
   * Enables the collection of information on subspaces during refinement
   */
  void enableStatsCollection();

  /**
   * Disables the collection of information on subspaces during refinement
   */
  void disableStatsCollection();

  /**
   * Returns information on all the subspaces that is sorted with respect to the iterations when
   * they have been added to the combigrid during refinement.
   *
   * NOTE: collectStats needs to be set in order to obtain useful information.
   *
   * @return the infos on the adaptive refinement steps
   */
  std::shared_ptr<LevelInfos> getInfoOnAddedLevels();
};
}  // namespace combigrid
}  // namespace sgpp
