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
#include <sgpp/combigrid/utils/DataVectorHashing.hpp>

#include <sgpp/combigrid/operation/multidim/LevelHelpers.hpp>  // TODO(holzmudd): remove with deprecated methods

#include <cmath>
#include <limits>
#include <memory>
#include <queue>
#include <string>
#include <unordered_set>
#include <vector>

namespace sgpp {
namespace combigrid {

/**
 * The CombigridEvaluator class evaluates a numerical method on different full grids and combines
 * them using the combination technique. The numerical method on a full grid is given by an
 * AbstractFullGridEvaluator.
 * The template parameter V determines whether this does single or multi evaluation, confer also the
 * description in algebraic/FloatArrayVector.hpp
 */
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

  /**
   * Evaluator that does the numerical evaluation on full grids.
   */
  std::shared_ptr<AbstractFullGridEvaluator<V>> multiEval;

  void initPartialDifferences() {
    partialDifferences.clear();
    for (size_t d = 0; d <= numDimensions; ++d) {
      partialDifferences.push_back(
          std::shared_ptr<AbstractMultiStorage<V>>(new TreeStorage<V>(numDimensions)));
    }
  }

 public:
  /**
   * Constructor.
   * @param numDimensions Dimensionality of the problem.
   * @param multiEval Evaluation method for full grids whose results are then combined into a single
   * value.
   */
  CombigridEvaluator(size_t numDimensions, std::shared_ptr<AbstractFullGridEvaluator<V>> multiEval)
      : sum(V::zero()),
        numDimensions(numDimensions),
        partialDifferences(),
        multiEval(multiEval),
        levelQueue(),
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
        levelQueue(),
        refinementStrategy(AdaptiveRefinementStrategy::minStrategy()) {
    initPartialDifferences();
  }

  /**
   * Performs an evaluation on the given level and adds the gained information to the current
   * numerical approximation value.
   * @return Returns true if the level was already there or (if it was computed and the difference
   * was not nan or +-inf).
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

  /**
   * @return a vector of tasks which can be precomputed in parallel to make the (serialized)
   * execution of addLevel() faster
   * @param level the level which one wants to compute
   * @param callback This callback is called (with already locked mutex) from inside one of the
   * returned tasks when all tasks for the given level are completed and the level can be added.
   */
  std::vector<ThreadPool::Task> getLevelTasks(MultiIndex const &level, ThreadPool::Task callback) {
    return multiEval->getLevelTasks(level, callback);
  }

  /**
   * Sets the mutex that is locked (if not nullptr) whenever problematic operations on data are
   * executed.
   */
  virtual void setMutex(std::shared_ptr<std::mutex> mutexPtr) { multiEval->setMutex(mutexPtr); }

  /**
   * Returns the storage with the Delta-values.
   */
  std::shared_ptr<AbstractMultiStorage<V>> differences() const {
    return partialDifferences[numDimensions];
  }

  /**
   * @returns true iff the given level has already been added.
   */
  bool containsLevel(MultiIndex const &level) {
    return partialDifferences[0]->containsIndex(level);
  }

  /**
   * Deprecated. Functionality has moved to LevelManager.
   */
  void setRefinementStrategy(AdaptiveRefinementStrategy strategy) { refinementStrategy = strategy; }

  /**
   * @return an estimate (upper bound, in the case of nested points normally exact) of the number of
   * new function evaluations (grid points) that have to be performed when adding this level.
   */
  size_t maxNewPoints(MultiIndex const &level) { return multiEval->maxNewPoints(level); }

  /**
   * @return the total number of grid points in a given level.
   */
  size_t numPoints(MultiIndex const &level) { return multiEval->numPoints(level); }

  /**
   * @return a vector with all grid points where the function has been evaluated (without
   * duplicates).
   */
  std::vector<base::DataVector> getAllGridPoints() {
    std::unordered_set<base::DataVector, DataVectorHash, DataVectorEqualTo> resultSet;

    auto it = partialDifferences[numDimensions]->getStoredDataIterator();

    while (it->isValid()) {
      auto fullGridPoints = multiEval->getGridPoints(it->getMultiIndex());

      for (auto &p : fullGridPoints) {
        resultSet.insert(p);
      }

      it->moveToNext();
    }

    std::vector<base::DataVector> result(resultSet.begin(), resultSet.end());

    return result;
  }

  /**
   * @return a DataMatrix with all grid points where the function has been evaluated in its columns
   * (without duplicates)
   */
  base::DataMatrix getGridPointMatrix() {
    auto vec = getAllGridPoints();
    base::DataMatrix m(vec[0].getSize(), vec.size());  // TODO(holzmudd): safety check

    for (size_t i = 0; i < vec.size(); ++i) {
      m.setColumn(i, vec[i]);
    }

    return m;
  }

  /**
   * @return the norm of the difference/Delta value (in the single-evaluation case this is the
   * absolute value) for the given level.
   */
  double getDifferenceNorm(MultiIndex const &level) {
    return partialDifferences[numDimensions]->get(level).norm();
  }

  /**
   * Deprecated. Functionality has moved to LevelManager.
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

  /**
   * Deprecated. Functionality has moved to LevelManager.
   */
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

  /**
   * @return the dimensionality of the problem.
   */
  size_t dim() const { return numDimensions; }

  /**
   * @return the numerical approximation value computed by the combination technique. No computation
   * is done here.
   */
  V getValue() const { return sum; }

  /**
   * Clears the already computed values. This method has to be called if a parameter changed etc.
   */
  void clear() {
    sum = V::zero();
    initPartialDifferences();
  }

  /**
   * @return a TreeStorage which contains an entry at an index i iff the level i has been added to
   * this CombigridEvaluator.
   */
  std::shared_ptr<TreeStorage<uint8_t>> getLevelStructure() {
    std::shared_ptr<TreeStorage<uint8_t>> storage(new TreeStorage<uint8_t>(numDimensions));

    auto it = partialDifferences[numDimensions]->getStoredDataIterator();

    while (it->isValid()) {
      storage->set(it->getMultiIndex(), 1);
      it->moveToNext();
    }

    return storage;
  }

  /**
   * @return the serialized version of getLevelStructure()
   */
  std::string getSerializedLevelStructure() {
    return TreeStorageSerializationStrategy<uint8_t>(numDimensions).serialize(getLevelStructure());
  }

  /**
   * Adds all levels for which an entry is contained in storage. "Inverse" operation to
   * getLevelStructure().
   */
  void addLevelsFromStructure(std::shared_ptr<TreeStorage<uint8_t>> storage) {
    auto it = storage->getStoredDataIterator();

    while (it->isValid()) {
      addLevel(it->getMultiIndex());
      it->moveToNext();
    }
  }

  /**
   * Equivalent to deserializing serializedStructure and then calling addLevelsFromStructure().
   * "Inverse" operation to getSerializedLevelStructure().
   */
  void addLevelsFromSerializedStructure(std::string serializedStructure) {
    addLevelsFromStructure(
        TreeStorageSerializationStrategy<uint8_t>(numDimensions).deserialize(serializedStructure));
  }

  /**
   * Deprecated. Functionality has moved to LevelManager.
   */
  void addLevelsAdaptive(size_t maxNumPoints) {
    CGLOG("addLevelsAdaptive(): start");
    initAdaption();

    CGLOG("addLevelsAdaptive(): adaption initialized");
    size_t currentPointBound = 0;

    while (!levelQueue.empty()) {
      CGLOG("addLevelsAdaptive(): next queue entry");
      QueueEntry entry = levelQueue.top();
      levelQueue.pop();

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
  // This whole part is deprecated. Functionality has moved to LevelManager.

  // data structures for adaptive refinement
  typedef std::priority_queue<QueueEntry, std::vector<QueueEntry>, QueueComparator> MultiIndexQueue;
  MultiIndexQueue levelQueue;
  AdaptiveRefinementStrategy refinementStrategy;

  void clearQueue() {
    MultiIndexQueue other;
    levelQueue.swap(other);  // unfortunately, there is no clear() method
  }

  void initAdaption() {
    clearQueue();

    auto it = partialDifferences[numDimensions]->getStoredDataIterator();

    while (it->isValid()) {
      tryAddSuccessorsToQueue(it->getMultiIndex());
      it->moveToNext();
    }

    if (levelQueue.empty()) {
      // no difference has yet been computed, so add start value
      MultiIndex startIndex(numDimensions, 0);
      levelQueue.push(QueueEntry(startIndex, 1.0, multiEval->maxNewPoints(startIndex)));
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
      levelQueue.push(nextEntry);
    }
  }
};
}  // namespace combigrid
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_MULTIDIM_COMBIGRIDEVALUATOR_HPP_ */
