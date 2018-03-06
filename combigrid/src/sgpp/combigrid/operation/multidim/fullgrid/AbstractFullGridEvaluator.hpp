// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/combigrid/definitions.hpp>
#include <sgpp/combigrid/grid/TensorGrid.hpp>
#include <sgpp/combigrid/grid/hierarchy/AbstractPointHierarchy.hpp>
#include <sgpp/combigrid/storage/AbstractCombigridStorage.hpp>

#include <memory>
#include <vector>

namespace sgpp {
namespace combigrid {

/**
 * Abstract class for doing an evaluation on a full grid, yielding a value of the template type V.
 * For this type, confer the comment in algebraic/FloatArrayVector.hpp
 * This class is used inside CombigridEvaluator to do evaluations on different grids.
 */
template <typename V>
class AbstractFullGridEvaluator {
 protected:
  /**
   * Provides access to the function values (stored or computed on demand)
   */
  std::shared_ptr<AbstractCombigridStorage> storage;
  // one per dimension
  std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies;

  /**
   * Pointer to a mutex that is locked when doing critical operations on data.
   * This is set to nullptr if the action is done in a single thread.
   */
  std::shared_ptr<std::recursive_mutex> mutexPtr;

 public:
  /**
   * Constructor.
   *
   * @param storage Storage that stores and provides the function values for each grid point.
   * @param pointHierarchies PointHierarchy objects for each dimension providing the points for each
   * level and information about their ordering.
   */
  AbstractFullGridEvaluator(std::shared_ptr<AbstractCombigridStorage> storage,
                            std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies)
      : storage(storage), pointHierarchies(pointHierarchies) {
    // TODO(holzmudd): check for dimension equality
  }

  virtual ~AbstractFullGridEvaluator() {}

  /**
   * Updates the current mutex. If the mutex is set to nullptr, no mutex locking is done. Otherwise,
   * the mutex is locked at critical actions.
   */
  virtual void setMutex(std::shared_ptr<std::recursive_mutex> mutexPtr) {
    this->mutexPtr = mutexPtr;
    storage->setMutex(mutexPtr);
  }

  /**
   * @return a vector of tasks which can be precomputed in parallel to make the (serialized)
   * execution of eval() faster
   * @param level the level which one wants to compute
   * @param callback This callback is called (with already locked mutex) from inside one of the
   * returned tasks when all tasks for the given level are completed and the level can be added.
   */
  virtual std::vector<ThreadPool::Task> getLevelTasks(MultiIndex const &level,
                                                      ThreadPool::Task callback) = 0;

  /**
   * Evaluates the function given through the storage for a certain level-multi-index (see class
   * description).
   */
  virtual V eval(MultiIndex const &level) = 0;

  /**
   * @return Returns the function value storage.
   */
  std::shared_ptr<AbstractCombigridStorage> getStorage() { return storage; }

  /**
   * Sets the parameters for the evaluators. Each dimension in which the evaluator does not need a
   * parameter is skipped.
   * So if only the evaluators at dimensions 1 and 3 need a parameter, params.size() should be 2 (or
   * at least 2)
   */
  virtual void setParameters(std::vector<V> const &params) = 0;

  /**
   * @return an estimate (upper bound, in the case of nested points normally exact) of the number of
   * new function evaluations (grid points) that have to be performed when evaluating this level.
   */
  virtual size_t maxNewPoints(MultiIndex const &level) {
    size_t result = 1;

    for (size_t d = 0; d < pointHierarchies.size(); ++d) {
      size_t currentLevel = level[d];
      size_t levelPoints = pointHierarchies[d]->getNumPoints(currentLevel);

      if (!pointHierarchies[d]->isNested() || currentLevel == 0) {
        result *= levelPoints;
      } else {
        result *= (levelPoints - pointHierarchies[d]->getNumPoints(currentLevel - 1));
      }
    }

    return result;
  }

  /**
   * @return Returns all grid points in the given level.
   */
  virtual std::vector<base::DataVector> getGridPoints(MultiIndex const &level) {
    size_t numDimensions = pointHierarchies.size();
    std::vector<base::DataVector> result;

    MultiIndex multiBounds(numDimensions);
    for (size_t d = 0; d < numDimensions; ++d) {
      multiBounds[d] = pointHierarchies[d]->getNumPoints(level[d]);
    }

    MultiIndexIterator it(multiBounds);

    while (it.isValid()) {
      auto index = it.getMultiIndex();

      base::DataVector vec(numDimensions);

      for (size_t d = 0; d < numDimensions; ++d) {
        vec[d] = pointHierarchies[d]->getPoint(level[d], index[d]);
      }

      result.push_back(vec);

      it.moveToNext();
    }

    return result;
  }

  /**
   * @return Returns the grid in the current level as a pointer to a TensorGrid object
   */
  virtual std::shared_ptr<TensorGrid> getTensorGrid(MultiIndex const &level,
                                                    std::vector<bool> orderingConfiguration) {
    size_t numDimensions = pointHierarchies.size();
    std::vector<base::DataVector> grids1D;

    for (size_t d = 0; d < numDimensions; ++d) {
      bool sorted = orderingConfiguration[d];
      grids1D.push_back(base::DataVector(pointHierarchies[d]->getPoints(level[d], sorted)));
    }

    return std::make_shared<TensorGrid>(grids1D, level);
  }

  /**
   * @return the total number of grid points in a given level.
   */
  virtual size_t numPoints(MultiIndex const &level) {
    size_t result = 1;

    for (size_t d = 0; d < pointHierarchies.size(); ++d) {
      size_t currentLevel = level[d];
      size_t levelPoints = pointHierarchies[d]->getNumPoints(currentLevel);
      result *= levelPoints;
    }

    return result;
  }
};

} /* namespace combigrid */
} /* namespace sgpp*/
