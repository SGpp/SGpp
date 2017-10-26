// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_MULTIDIM_FULLGRID_ABSTRACTFULLGRIDEVALUATIONSTRATEGY_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_MULTIDIM_FULLGRID_ABSTRACTFULLGRIDEVALUATIONSTRATEGY_HPP_

#include <sgpp/combigrid/GeneralFunction.hpp>
#include <sgpp/combigrid/grid/TensorGrid.hpp>
#include <sgpp/combigrid/grid/hierarchy/AbstractPointHierarchy.hpp>
#include <sgpp/combigrid/operation/multidim/fullgrid/AbstractFullGridEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/AbstractLinearEvaluator.hpp>
#include <sgpp/combigrid/storage/AbstractCombigridStorage.hpp>
#include <sgpp/combigrid/storage/tree/TreeStorage.hpp>
#include <sgpp/combigrid/threading/PtrGuard.hpp>

#include <vector>

namespace sgpp {
namespace combigrid {
/**
 * currently two evaluation strategies are suppoted:
 * linear: uses the eval function of AbstractFullGridLinearEvaluator, i.e. evaluation in every
 *		   point individually
 * grid_based: uses a gridFunction to evaluate in all grid points at once
 *
 * You decide which one to use simply by calling the constructor with or without a gridFunction
*/
enum class Strategy { linear, grid_based };

typedef GeneralFunction<std::shared_ptr<TreeStorage<double>>, std::shared_ptr<TensorGrid>>
    GridFunction;

template <typename V>
class AbstractFullGridEvaluationStrategy : public AbstractFullGridEvaluator<V> {
 public:
  AbstractFullGridEvaluationStrategy(
      std::shared_ptr<AbstractCombigridStorage> storage,
      std::vector<std::shared_ptr<AbstractLinearEvaluator<V>>> evaluatorPrototypes,
      std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies)
      : AbstractFullGridEvaluator<V>(storage, pointHierarchies),
        storage(storage),
        evaluatorPrototypes(evaluatorPrototypes),
        pointHierarchies(pointHierarchies),
        gridFunction(nullptr),
        strategy(Strategy::linear){};
  AbstractFullGridEvaluationStrategy(
      std::shared_ptr<AbstractCombigridStorage> storage,
      std::vector<std::shared_ptr<AbstractLinearEvaluator<V>>> evaluatorPrototypes,
      std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
      GridFunction gridFunction)
      : AbstractFullGridEvaluator<V>(storage, pointHierarchies),
        storage(storage),
        evaluatorPrototypes(evaluatorPrototypes),
        pointHierarchies(pointHierarchies),
        gridFunction(gridFunction),
        strategy(Strategy::grid_based),
        precomputedLevels(std::make_shared<TreeStorage<uint8_t>>(pointHierarchies.size())){};

  ~AbstractFullGridEvaluationStrategy(){};

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
  /**
   * @return a vector of tasks which can be precomputed in parallel to make the (serialized)
   * execution of eval() faster
   * @param level the level which one wants to compute
   * @param callback This callback is called (with already locked mutex) from inside one of the
   * returned tasks when all tasks for the given level are completed and the level can be added.
   */
  std::vector<ThreadPool::Task> getLevelTasksLinear(MultiIndex const &level,
                                                    ThreadPool::Task callback) {
    size_t numDimensions = this->evaluatorPrototypes.size();
    MultiIndex multiBounds(numDimensions);

    for (size_t d = 0; d < numDimensions; ++d) {
      multiBounds[d] = this->pointHierarchies[d]->getNumPoints(level[d]);
    }

    MultiIndexIterator it(multiBounds);
    std::vector<bool> orderingConfiguration(numDimensions,
                                            false);  // The request does not need ordered points
    auto funcIter = this->storage->getGuidedIterator(level, it, orderingConfiguration);

    std::vector<std::function<double()>> computationTasks;
    std::vector<MultiIndex> multiIndices;

    while (funcIter->isValid()) {
      if (!funcIter->computationRequested()) {
        computationTasks.push_back(funcIter->requestComputationTask());
        multiIndices.push_back(funcIter->getMultiIndex());
      }
      funcIter->moveToNext();
    }

    std::vector<ThreadPool::Task> tasks;

    if (computationTasks.empty()) {
      callback();
    } else {
      // make it a pointer so that it does not get deleted before all tasks are completed
      auto counter = std::make_shared<size_t>(computationTasks.size());

      for (size_t i = 0; i < computationTasks.size(); ++i) {
        auto compTask = computationTasks[i];
        auto index = multiIndices[i];

        tasks.push_back(ThreadPool::Task([compTask, index, counter, callback, this, level]() {
          auto result = compTask();

          CGLOG_SURROUND(PtrGuard guard(this->mutexPtr));
          this->storage->set(level, index, result);
          --(*counter);
          if (*counter == 0) {
            callback();
          }
          CGLOG("leave guard(this->mutexPtr) in FGEval");
        }));
      }
    }

    return tasks;
  }

  std::shared_ptr<TensorGrid> getTensorGrid2(MultiIndex const &level) {
    std::vector<bool> orderingConfiguration(this->evaluatorPrototypes.size());
    for (size_t d = 0; d < this->evaluatorPrototypes.size(); ++d) {
      orderingConfiguration[d] = this->evaluatorPrototypes[d]->needsOrderedPoints();
    }
    return this->getTensorGrid(level, orderingConfiguration);
  }

  /**
     * @return a vector of tasks which can be precomputed in parallel to make the (serialized)
     * execution of eval() faster. This class only returns one task in the vector.
     * @param level the level which one wants to compute
     * @param callback This callback is called (with already locked mutex) from inside one of the
     * returned tasks when all tasks for the given level are completed and the level can be added.
     */
  std::vector<ThreadPool::Task> getLevelTasksGridBased(MultiIndex const &level,
                                                       ThreadPool::Task callback) {
    auto grid = getTensorGrid2(level);

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

  std::vector<ThreadPool::Task> getLevelTasks(MultiIndex const &level,
                                              ThreadPool::Task callback) override {
    if (strategy == Strategy::linear)
      return getLevelTasksLinear(level, callback);
    else  // if (strategy == Strategy::grid_based)
      return getLevelTasksGridBased(level, callback);
  }

  //  V eval(MultiIndex const &level) override {
  //    if ((!precomputedLevels->containsIndex(level)) && (strategy == Strategy::grid_based)) {
  //        addResults(level, gridFunction(getTensorGrid2(level)));
  //      precomputedLevels->set(level, 1);
  //    }
  //    // call the base eval
  //    return combinationEvaluator->eval(level);
  //}

  //  void setParameters(std::vector<V> const &params) {
  //  combinationEvaluator->setParameters(params);
  //  }

 protected:
  std::shared_ptr<AbstractCombigridStorage> storage;
  std::vector<std::shared_ptr<AbstractLinearEvaluator<V>>> evaluatorPrototypes;
  std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies;
  GridFunction gridFunction;
  Strategy strategy;
  // We cannot use bool because that would involve vector<bool> problems...
  std::shared_ptr<TreeStorage<uint8_t>> precomputedLevels;
};

} /* namespace combigrid */
} /* namespace sgpp */

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_MULTIDIM_FULLGRID_ABSTRACTFULLGRIDEVALUATIONSTRATEGY_HPP_ */
