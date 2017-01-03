// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

namespace sgpp {
namespace combigrid {

#include <sgpp/combigrid/common/MultiIndexIterator.hpp>
#include <sgpp/combigrid/definitions.hpp>
#include <sgpp/combigrid/grid/hierarchy/AbstractPointHierarchy.hpp>
#include <sgpp/combigrid/operation/multidim/AbstractFullGridLinearEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/AbstractLinearEvaluator.hpp>
#include <sgpp/combigrid/storage/AbstractCombigridStorage.hpp>
#include <sgpp/combigrid/threading/PtrGuard.hpp>
#include <sgpp/combigrid/threading/ThreadPool.hpp>

#include <vector>

template <typename V>
class FullGridLinearCallbackEvaluator : public AbstractFullGridLinearEvaluator<V> {
 public:
  /**
   * Constructor.
   *
   * @param storage Storage that stores and provides the function values for each grid point.
   * @param evaluatorPrototypes prototype objects for the evaluators that are cloned to get an
   * evaluator for each dimension and each level.
   * @param pointHierarchies PointHierarchy objects for each dimension providing the points for each
   * level and information about their ordering.
   */
  FullGridLinearCallbackEvaluator(
      std::shared_ptr<AbstractCombigridStorage> storage,
      std::vector<std::shared_ptr<AbstractLinearEvaluator<V>>> evaluatorPrototypes,
      std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies)
      : AbstractFullGridTensorEvaluator(storage, evaluatorPrototypes, pointHierarchies) {}

  virtual ~FullGridLinearCallbackEvaluator() {}

  /**
   * @return a vector of tasks which can be precomputed in parallel to make the (serialized)
   * execution of eval() faster
   * @param level the level which one wants to compute
   * @param callback This callback is called (with already locked mutex) from inside one of the
   * returned tasks when all tasks for the given level are completed and the level can be added.
   */
  std::vector<ThreadPool::Task> getLevelTasks(MultiIndex const &level, ThreadPool::Task callback) {
    size_t numDimensions = evaluators.size();
    MultiIndex multiBounds(numDimensions);

    for (size_t d = 0; d < numDimensions; ++d) {
      multiBounds[d] = pointHierarchies[d]->getNumPoints(level[d]);
    }

    MultiIndexIterator it(multiBounds);
    std::vector<bool> orderingConfiguration(numDimensions,
                                            false);  // The request does not need ordered points
    auto funcIter = storage->getGuidedIterator(level, it, orderingConfiguration);

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

        tasks.push_back([compTask, index, counter, callback, this, level]() {
          auto result = compTask();

          CGLOG_SURROUND(PtrGuard guard(this->mutexPtr));
          this->storage->set(level, index, result);
          --(*counter);
          if (*counter == 0) {
            callback();
          }
          CGLOG("leave guard(this->mutexPtr) in FGEval");
        });
      }
    }

    return tasks;
  }
};

} /* namespace combigrid */
} /* namespace sgpp */
