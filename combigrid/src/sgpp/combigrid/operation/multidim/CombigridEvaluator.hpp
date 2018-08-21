// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_MULTIDIM_COMBIGRIDEVALUATOR_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_MULTIDIM_COMBIGRIDEVALUATOR_HPP_

#include <sgpp/combigrid/common/BoundedSumMultiIndexIterator.hpp>
#include <sgpp/combigrid/definitions.hpp>
#include <sgpp/combigrid/operation/multidim/AbstractLevelEvaluator.hpp>
#include <sgpp/combigrid/operation/multidim/AdaptiveRefinementStrategy.hpp>
#include <sgpp/combigrid/operation/multidim/fullgrid/AbstractFullGridEvaluationStrategy.hpp>
#include <sgpp/combigrid/serialization/TreeStorageSerializationStrategy.hpp>
#include <sgpp/combigrid/storage/AbstractMultiStorage.hpp>
#include <sgpp/combigrid/storage/tree/TreeStorage.hpp>
#include <sgpp/combigrid/utils/DataVectorHashing.hpp>
#include <sgpp/combigrid/algebraic/NormStrategy.hpp>

#include <cmath>
#include <iostream>
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

  /**
   * strategy that computes the norm of the full grid evaluation results
   */
  std::shared_ptr<NormStrategy<V>> normStrategy;

  /**
   * Upper bound for the number of points (function evaluations) used.
   */
  size_t upperPointBound = 0;

  void initPartialDifferences() {
    partialDifferences.clear();
    for (size_t d = 0; d <= numDimensions; ++d) {
      partialDifferences.push_back(std::make_shared<TreeStorage<V>>(numDimensions));
    }
  }

 public:
  /**
   * Constructor.
   * @param numDimensions Dimensionality of the problem.
   * @param multiEval Evaluation method for full grids whose results are then combined into a single
   * value.
   * @param normStrategy defines how the norm of differences is computed.
   */
  CombigridEvaluator(size_t numDimensions, std::shared_ptr<AbstractFullGridEvaluator<V>> multiEval,
                     std::shared_ptr<NormStrategy<V>> normStrategy = nullptr)
      : sum(V::zero()),
        numDimensions(numDimensions),
        partialDifferences(),
        multiEval(multiEval),
        normStrategy(normStrategy) {
    if (this->normStrategy == nullptr) {
      this->normStrategy = std::make_shared<NormStrategy<V>>();
    }
    initPartialDifferences();
  }

  // ToDo (rehmemk) bei Umstellung von AbstractFullGridLinearEvaluator auf Summation- und
  // EvalStrategy das folgende einfach auskommentiert. Prüfen und gegebenenfalls wieder einführen!

  /**
   * For some reason, SWIG cannot convert the shared_ptr into the more abstract type, so we need
   * this 'duplicate'
   */
  //  CombigridEvaluator(size_t numDimensions,
  //                     std::shared_ptr<FullGridLinearCallbackEvaluator<V>> multiEval)
  //      : sum(V::zero()), numDimensions(numDimensions), partialDifferences(), multiEval(multiEval)
  //      {
  //    initPartialDifferences();
  //  }

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
    upperPointBound += maxNewPoints(level);
    CGLOG("addLevel(): add previous levels");
    // ensure that preceding indices are already computed
    for (size_t d = 0; d < numDimensions; ++d) {
      if (level[d] > 0) {
        MultiIndex l = level;
        --l[d];

        // should not affect performance because nothing is computed if the storage
        // already contains a value
        bool success = addLevel(l);
        if (!success) {
          return false;
        }
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
    double norm = normStrategy->norm(value);
    if (std::isnan(norm) || std::isinf(norm)) {
      return false;
    } else {
      sum.add(value);
      return true;
    }
  }

  /**
   * @return An upper bound for the number of points (function evaluations) used for the current
   * computation. This bound is exact if nesting is used or if otherwise each grid point only occurs
   * in exactly one level.
   */
  size_t getUpperPointBound() const { return upperPointBound; }

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
  virtual void setMutex(std::shared_ptr<std::recursive_mutex> mutexPtr) {
    multiEval->setMutex(mutexPtr);
  }

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
   * @return an estimate (upper bound, in the case of nested points normally exact) of the number of
   * new function evaluations (grid points) that have to be performed when adding this level.
   */
  size_t maxNewPoints(MultiIndex const &level) { return multiEval->maxNewPoints(level); }

  size_t maxNumPointsForRegular(size_t q) {
    auto it = std::make_shared<BoundedSumMultiIndexIterator>(numDimensions, q);
    size_t sum = 0;
    for (; it->isValid(); it->moveToNext()) {
      sum += maxNewPoints(it->value());
    }
    return sum;
  }

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
    return normStrategy->norm(partialDifferences[numDimensions]->get(level));
  }

  /**
   * @return the dimensionality of the problem.
   */
  size_t numDims() const { return numDimensions; }

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
    upperPointBound = 0;
    initPartialDifferences();
  }

  /**
   * @return a TreeStorage which contains an entry at an index i iff the level i has been added to
   * this CombigridEvaluator.
   */
  std::shared_ptr<TreeStorage<uint8_t>> getLevelStructure() {
    std::shared_ptr<TreeStorage<uint8_t>> storage =
        std::make_shared<TreeStorage<uint8_t>>(numDimensions);

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
};
}  // namespace combigrid
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_MULTIDIM_COMBIGRIDEVALUATOR_HPP_ */
