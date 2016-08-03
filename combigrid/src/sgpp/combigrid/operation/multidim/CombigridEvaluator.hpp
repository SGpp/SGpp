// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_MULTIDIM_COMBIGRIDEVALUATOR_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_MULTIDIM_COMBIGRIDEVALUATOR_HPP_

#include <sgpp/combigrid/common/BoundedSumMultiIndexIterator.hpp>
#include <sgpp/combigrid/definitions.hpp>
#include <sgpp/combigrid/operation/multidim/AbstractFullGridEvaluator.hpp>
#include <sgpp/combigrid/operation/multidim/AbstractLevelEvaluator.hpp>
#include <sgpp/combigrid/operation/multidim/AdaptiveRefinementStrategy.hpp>
#include <sgpp/combigrid/operation/multidim/FullGridTensorEvaluator.hpp>
#include <sgpp/combigrid/serialization/TreeStorageSerializationStrategy.hpp>
#include <sgpp/combigrid/storage/AbstractMultiStorage.hpp>
#include <sgpp/combigrid/storage/tree/TreeStorage.hpp>

#include <sgpp/combigrid/operation/multidim/LevelHelpers.hpp>  // TODO(holzmuedd): remove

#include <cmath>
#include <limits>
#include <memory>
#include <queue>
#include <string>
#include <unordered_set>
#include <vector>

namespace sgpp {
namespace combigrid {

template <typename V>
class CombigridEvaluator : public AbstractLevelEvaluator {
  V sum;
  size_t numDimensions;

  /**
   * All partial differences are stored in these storages.
   * I. e. storages[0] contains the values for the full grids
   * and storages[numDimensions] contains the normal differences.
   * The storages in between contain partial differences.
   */
  std::vector<std::shared_ptr<AbstractMultiStorage<V>>> partialDifferences;
  std::shared_ptr<AbstractFullGridEvaluator<V>> multiEval;

  void initPartialDifferences() {
    partialDifferences.clear();
    for (size_t d = 0; d <= numDimensions; ++d) {
      partialDifferences.push_back(
          std::shared_ptr<AbstractMultiStorage<V>>(new TreeStorage<V>(numDimensions)));
    }
  }

 public:
  CombigridEvaluator(size_t numDimensions, std::shared_ptr<AbstractFullGridEvaluator<V>> multiEval)
      : sum(V::zero()),
        numDimensions(numDimensions),
        partialDifferences(),
        multiEval(multiEval),
        queue(),
        refinementStrategy(AdaptiveRefinementStrategy::minStrategy()) {
    initPartialDifferences();
  }

  /**
   * For some reason, SWIG cannot convert the shared_ptr into the more abstract type, so we need
   * this 'duplicate'
   */
  CombigridEvaluator(size_t numDimensions, std::shared_ptr<FullGridTensorEvaluator<V>> multiEval)
      : sum(V::zero()),
        numDimensions(numDimensions),
        partialDifferences(),
        multiEval(multiEval),
        queue(),
        refinementStrategy(AdaptiveRefinementStrategy::minStrategy()) {
    initPartialDifferences();
  }

  /**
   * @return Returns true if the level was already there or if it was computed and the difference
   * was not nan or +-inf.
   */
  bool addLevel(MultiIndex const &level) {
    CGLOG("addLevel(): start");
    if (containsLevel(level)) {
      return true;
    }

    CGLOG("addLevel(): add previous levels");
    // ensure that preceding indices are already computed
    for (size_t d = 0; d < numDimensions; ++d) {
      if (level[d] > 0) {
        MultiIndex l = level;
        --l[d];
        addLevel(l);  // should not affect performance because nothing is computed if the storage
                      // already contains a value
      }
    }

    CGLOG("addLevel(): eval level");
    V value = multiEval->eval(level);

    CGLOG("addLevel(): evaluate partial differences");
    partialDifferences[0]->set(level, value);

    for (size_t d = 0; d < numDimensions; ++d) {
      if (level[d] == 0) {
        partialDifferences[d + 1]->set(level, value);
      } else {
        MultiIndex predecessor = level;
        --predecessor[d];
        value.sub(partialDifferences[d]->get(predecessor));
        partialDifferences[d + 1]->set(level, value);
      }
    }

    double norm = value.norm();

    if (std::isnan(norm) || std::isinf(norm)) {
      return false;
    } else {
      sum.add(value);
      return true;
    }
  }

  std::vector<ThreadPool::Task> getLevelTasks(MultiIndex const &level, ThreadPool::Task callback) {
    return multiEval->getLevelTasks(level, callback);
  }

  virtual void setMutex(std::shared_ptr<std::mutex> mutexPtr) { multiEval->setMutex(mutexPtr); }

  std::shared_ptr<AbstractMultiStorage<V>> differences() const {
    return partialDifferences[numDimensions];
  }

  bool containsLevel(MultiIndex const &level) {
    return partialDifferences[0]->containsIndex(level);
  }

  void setRefinementStrategy(AdaptiveRefinementStrategy strategy) { refinementStrategy = strategy; }

  size_t maxNewPoints(MultiIndex const &level) { return multiEval->maxNewPoints(level); }

  size_t numPoints(MultiIndex const &level) { return multiEval->numPoints(level); }

  double getDifferenceNorm(MultiIndex const &level) {
    return partialDifferences[numDimensions]->get(level).norm();
  }

  /**
   * @param q: Maximum 1-norm of the level-multi-index, where the levels start from 0 (not from 1 as
   * in most papers).
   * If you have a norm w with levels starting from 1, simply use q = w - dim().
   */
  void addRegularLevels(size_t q) {
    BoundedSumMultiIndexIterator iterator(numDimensions, q);

    while (iterator.isValid()) {
      addLevel(iterator.value());
      iterator.moveToNext();
    }
  }

  void addRegularLevelsByNumPoints(size_t maxNumPoints) {
    size_t numPoints = 0;

    size_t q = 0;

    while (true) {
      BoundedSumMultiIndexIterator iterator(numDimensions, q);

      while (iterator.isValid()) {
        MultiIndex nextLevel = iterator.value();

        if (partialDifferences[0]->containsIndex(nextLevel)) {
          iterator.moveToNext();
          continue;
        }

        size_t maxNewPoints = multiEval->maxNewPoints(nextLevel);
        numPoints += maxNewPoints;

        if (numPoints > maxNumPoints) {
          return;
        }

        addLevel(nextLevel);
        iterator.moveToNext();
      }

      ++q;
    }
  }

  size_t dim() const { return numDimensions; }

  V getValue() const { return sum; }

  void clear() {
    sum = V::zero();
    initPartialDifferences();
  }

  std::shared_ptr<TreeStorage<uint8_t>> getLevelStructure() {
    std::shared_ptr<TreeStorage<uint8_t>> storage(new TreeStorage<uint8_t>(numDimensions));

    auto it = partialDifferences[numDimensions]->getStoredDataIterator();

    while (it->isValid()) {
      storage->set(it->getMultiIndex(), 1);
      it->moveToNext();
    }

    return storage;
  }

  std::string getSerializedLevelStructure() {
    return TreeStorageSerializationStrategy<uint8_t>(numDimensions).serialize(getLevelStructure());
  }

  void addLevelsFromStructure(std::shared_ptr<TreeStorage<uint8_t>> storage) {
    auto it = storage->getStoredDataIterator();

    while (it->isValid()) {
      addLevel(it->getMultiIndex());
      it->moveToNext();
    }
  }

  void addLevelsFromSerializedStructure(std::string serializedStructure) {
    addLevelsFromStructure(
        TreeStorageSerializationStrategy<uint8_t>(numDimensions).deserialize(serializedStructure));
  }

  void addLevelsAdaptive(size_t maxNumPoints) {
    CGLOG("addLevelsAdaptive(): start");
    initAdaption();

    CGLOG("addLevelsAdaptive(): adaption initialized");
    size_t currentPointBound = 0;

    while (!queue.empty()) {
      CGLOG("addLevelsAdaptive(): next queue entry");
      QueueEntry entry = queue.top();
      queue.pop();

      currentPointBound += entry.maxNewPoints;

      if (currentPointBound > maxNumPoints) {
        break;
      }

      CGLOG("addLevelsAdaptive(): add next level");
      bool validResult = addLevel(entry.level);

      CGLOG("addLevelsAdaptive(): added next level");
      if (validResult) {
        tryAddSuccessorsToQueue(entry.level);
      }
    }
  }

 private:
  // data structures for adaptive refinement
  typedef std::priority_queue<QueueEntry, std::vector<QueueEntry>, QueueComparator> MultiIndexQueue;
  MultiIndexQueue queue;
  AdaptiveRefinementStrategy refinementStrategy;

  void clearQueue() {
    MultiIndexQueue other;
    queue.swap(other);  // unfortunately, there is no clear() method
  }

  void initAdaption() {
    clearQueue();

    auto it = partialDifferences[numDimensions]->getStoredDataIterator();

    while (it->isValid()) {
      tryAddSuccessorsToQueue(it->getMultiIndex());
      it->moveToNext();
    }

    if (queue.empty()) {
      // no difference has yet been computed, so add start value
      MultiIndex startIndex(numDimensions, 0);
      queue.push(QueueEntry(startIndex, 1.0, multiEval->maxNewPoints(startIndex)));
    }
  }

  void tryAddSuccessorsToQueue(MultiIndex const &index) {
    for (size_t d = 0; d < numDimensions; ++d) {
      MultiIndex nextIndex = index;
      ++nextIndex[d];

      tryAddIndexToQueue(nextIndex);
    }
  }

  void tryAddIndexToQueue(MultiIndex const &nextIndex) {
    if (partialDifferences[numDimensions]->containsIndex(nextIndex)) {
      return;
    }

    // calculate estimated value...
    std::vector<double> predecessorNorms;

    bool hasAllPredecessors = true;

    for (size_t prevDim = 0; prevDim < numDimensions; ++prevDim) {
      if (nextIndex[prevDim] > 0) {
        MultiIndex prevIndex = nextIndex;
        --prevIndex[prevDim];

        if (partialDifferences[numDimensions]->containsIndex(prevIndex)) {
          double norm = partialDifferences[numDimensions]->get(prevIndex).norm();
          predecessorNorms.push_back(norm);
        } else {
          hasAllPredecessors = false;  // add index later, not all predecessors are available
          break;
        }
      }
    }

    if (hasAllPredecessors) {
      size_t maxNewPoints = multiEval->maxNewPoints(nextIndex);
      double priority = refinementStrategy.computePriority(predecessorNorms, maxNewPoints);

      QueueEntry nextEntry(nextIndex, priority, maxNewPoints);
      queue.push(nextEntry);
    }
  }
};
}  // namespace combigrid
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_MULTIDIM_COMBIGRIDEVALUATOR_HPP_ */
