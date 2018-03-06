// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/combigrid/GeneralFunction.hpp>
#include <sgpp/combigrid/common/MultiIndexIterator.hpp>
#include <sgpp/combigrid/definitions.hpp>
#include <sgpp/combigrid/grid/TensorGrid.hpp>
#include <sgpp/combigrid/grid/hierarchy/AbstractPointHierarchy.hpp>
#include <sgpp/combigrid/operation/multidim/fullgrid/AbstractFullGridEvaluationStrategy.hpp>
#include <sgpp/combigrid/operation/onedim/AbstractLinearEvaluator.hpp>
#include <sgpp/combigrid/storage/AbstractCombigridStorage.hpp>
#include <sgpp/combigrid/storage/tree/TreeStorage.hpp>
#include <sgpp/combigrid/threading/PtrGuard.hpp>
#include <sgpp/combigrid/threading/ThreadPool.hpp>

#include <iostream>
#include <vector>

namespace sgpp {
namespace combigrid {

/**
 * Implementation of the AbstractFullGridLinearEvaluator class using a callback function that
 * If you want to be able to use different function values at the same point in different levels
 * (for example because you are implementing a PDE solver), set exploitNesting to false in the
 * constructor of CombigridTreeStorage.
 */
template <typename V>
class FullGridGridBasedEvaluator : public AbstractFullGridEvaluationStrategy<V> {
 protected:
  GridFunction gridFunction;

  // We cannot use bool because that would involve vector<bool> problems...
  std::shared_ptr<TreeStorage<uint8_t>> precomputedLevels;

  void addResults(MultiIndex const &level, std::shared_ptr<TreeStorage<double>> results) {
    std::vector<bool> orderingConfiguration(this->evaluatorPrototypes.size());
    for (size_t d = 0; d < this->evaluatorPrototypes.size(); ++d) {
      orderingConfiguration[d] = this->evaluatorPrototypes[d]->needsOrderedPoints();
    }
    MultiIndex multiBounds(this->pointHierarchies.size());
    for (size_t d = 0; d < multiBounds.size(); ++d) {
      multiBounds[d] = this->pointHierarchies[d]->getNumPoints(level[d]);
    }
    MultiIndexIterator multiIt(multiBounds);

    auto storageIt = this->storage->getGuidedIterator(level, multiIt, orderingConfiguration);

    while (storageIt->isValid()) {
      storageIt->value() = results->get(multiIt.getMultiIndex());
      storageIt->moveToNext();
    }
  }

 public:
  /**
   * Constructor.
   *
   * @param storage Storage that stores and provides the function values for each grid point.
   * @param evaluatorPrototypes prototype objects for the evaluators that are cloned to get an
   * evaluator for each dimension and each level.
   * @param pointHierarchies PointHierarchy objects for each dimension providing the points for
   * each level and information about their ordering.
   * @param gridFunction callback function that is called with a grid as parameters and should
   * @param summationStrategyType strategy to gather the results of the univariate evaluators on
   * each anisotropic full grid
   *
   * return a TreeStorage that contains the values at these grid points
   */
  FullGridGridBasedEvaluator(
      std::shared_ptr<AbstractCombigridStorage> storage,
      std::vector<std::shared_ptr<AbstractLinearEvaluator<V>>> evaluatorPrototypes,
      std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
      GridFunction gridFunction,
      FullGridSummationStrategyType summationStrategyType = FullGridSummationStrategyType::LINEAR)
      : AbstractFullGridEvaluationStrategy<V>(storage, evaluatorPrototypes, pointHierarchies,
                                              summationStrategyType),
        gridFunction(gridFunction),
        precomputedLevels(std::make_shared<TreeStorage<uint8_t>>(pointHierarchies.size())) {}

  virtual ~FullGridGridBasedEvaluator() {}

  /**
     * @return a vector of tasks which can be precomputed in parallel to make the (serialized)
     * execution of eval() faster. This class only returns one task in the vector.
     * @param level the level which one wants to compute
     * @param callback This callback is called (with already locked mutex) from inside one of the
     * returned tasks when all tasks for the given level are completed and the level can be added.
     */
  std::vector<ThreadPool::Task> getLevelTasks(MultiIndex const &level, ThreadPool::Task callback) {
    size_t numDimensions = this->pointHierarchies.size();
    std::vector<bool> orderingConfiguration(numDimensions);

    for (size_t d = 0; d < numDimensions; ++d) {
      bool needsOrdered = this->evaluatorPrototypes[d]->needsOrderedPoints();
      orderingConfiguration[d] = needsOrdered;
    }
    auto grid = this->getTensorGrid(level, orderingConfiguration);
    std::vector<ThreadPool::Task> tasks;

    tasks.push_back(ThreadPool::Task([grid, level, this, callback]() {
      auto results = gridFunction(grid);

      // now we need locking
      PtrGuard guard(this->mutexPtr);
      addResults(level, results);
      precomputedLevels->set(level, 1);
      callback();
    }));

    return tasks;
  }

  virtual V eval(MultiIndex const &level) {
    if (!precomputedLevels->containsIndex(level)) {
      size_t numDimensions = this->pointHierarchies.size();
      std::vector<bool> orderingConfiguration(numDimensions);

      for (size_t d = 0; d < numDimensions; ++d) {
        bool needsOrdered = this->evaluatorPrototypes[d]->needsOrderedPoints();
        orderingConfiguration[d] = needsOrdered;
      }
      addResults(level, gridFunction(this->getTensorGrid(level, orderingConfiguration)));

      precomputedLevels->set(level, 1);
    }

    // call the base eval
    return this->summationStrategy->eval(level);
  }
};

} /* namespace combigrid */
} /* namespace sgpp */
