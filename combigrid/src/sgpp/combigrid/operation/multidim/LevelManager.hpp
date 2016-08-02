// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_MULTIDIM_LEVELMANAGER_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_MULTIDIM_LEVELMANAGER_HPP_

#include <sgpp/globaldef.hpp>

#include <sgpp/combigrid/common/BoundedSumMultiIndexIterator.hpp>
#include <sgpp/combigrid/definitions.hpp>
#include <sgpp/combigrid/operation/multidim/AbstractLevelEvaluator.hpp>  // TODO(holzmudd): remove
#include <sgpp/combigrid/operation/multidim/AdaptiveRefinementStrategy.hpp>
#include <sgpp/combigrid/operation/multidim/LevelHelpers.hpp>
#include <sgpp/combigrid/serialization/TreeStorageSerializationStrategy.hpp>
#include <sgpp/combigrid/storage/AbstractMultiStorage.hpp>
#include <sgpp/combigrid/storage/tree/TreeStorage.hpp>

#include <cmath>
#include <limits>
#include <memory>
#include <mutex>
#include <queue>
#include <string>
#include <unordered_set>
#include <vector>

namespace sgpp {
namespace combigrid {

// TODO(holzmudd): put functions in cpp file

class LevelManager {
 protected:
  // data structures for adaptive refinement
  MultiIndexQueue queue;
  std::shared_ptr<TreeStorage<std::shared_ptr<LevelInfo>>> levelData;
  size_t numDimensions;
  std::shared_ptr<AbstractLevelEvaluator> combiEval;
  std::shared_ptr<std::mutex> managerMutex;

  /**
   * By implementing this method in a derived class, the adaption can be customized.
   */
  virtual double computePriority(MultiIndex const &level) = 0;

  /**
   * Initializes the data structures for adaptive level generation
   */
  virtual void initAdaption();

  virtual void tryAddSuccessors(MultiIndex const &level);

  virtual void tryAddLevel(MultiIndex const &level);

  virtual void addToQueue(MultiIndex const &level, std::shared_ptr<LevelInfo> levelInfo);

  virtual std::vector<MultiIndex> getPredecessors(MultiIndex const &level);

  virtual std::vector<MultiIndex> getSuccessors(MultiIndex const &level);

  virtual void beforeComputation(MultiIndex const &level);

  virtual void afterComputation(MultiIndex const &level);

  virtual void predecessorsCompleted(MultiIndex const &level);

  virtual void updatePriority(MultiIndex const &level, std::shared_ptr<LevelInfo> levelInfo);

  /**
   * @param q: Maximum 1-norm of the level-multi-index, where the levels start from 0 (not from 1 as
   * in most papers).
   * If you have a norm w with levels starting from 1, simply use q = w - dim().
   */
  std::vector<MultiIndex> getRegularLevels(size_t q);

  std::vector<MultiIndex> getRegularLevelsByNumPoints(size_t maxNumPoints);

  void precomputeLevelsParallel(std::vector<MultiIndex> const &levels, size_t numThreads);

  void addLevels(std::vector<MultiIndex> const &levels);

 public:
  LevelManager(std::shared_ptr<AbstractLevelEvaluator> levelEvaluator);

  virtual ~LevelManager();

  /**
   * @param q: Maximum 1-norm of the level-multi-index, where the levels start from 0 (not from 1 as
   * in most papers).
   * If you have a norm w with levels starting from 1, simply use q = w - dim().
   */
  void addRegularLevels(size_t q);

  void addRegularLevelsByNumPoints(size_t maxNumPoints);

  /**
   * @param q: Maximum 1-norm of the level-multi-index, where the levels start from 0 (not from 1 as
   * in most papers).
   * If you have a norm w with levels starting from 1, simply use q = w - dim().
   */
  void addRegularLevelsParallel(size_t q, size_t numThreads);

  void addRegularLevelsByNumPointsParallel(size_t maxNumPoints, size_t numThreads);

  virtual size_t dim() const;

  std::shared_ptr<TreeStorage<uint8_t>> getLevelStructure() const;

  std::string getSerializedLevelStructure() const;

  void addLevelsFromStructure(std::shared_ptr<TreeStorage<uint8_t>> storage);

  void addLevelsFromSerializedStructure(std::string serializedStructure);

  virtual void addLevelsAdaptive(size_t maxNumPoints);

  virtual void addLevelsAdaptiveParallel(size_t maxNumPoints, size_t numThreads);
};
}  // namespace combigrid
}  // namespace sgpp

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_MULTIDIM_LEVELMANAGER_HPP_ */
