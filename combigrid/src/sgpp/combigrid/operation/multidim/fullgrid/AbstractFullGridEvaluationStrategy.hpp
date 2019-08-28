// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/combigrid/GeneralFunction.hpp>
#include <sgpp/combigrid/grid/TensorGrid.hpp>
#include <sgpp/combigrid/grid/hierarchy/AbstractPointHierarchy.hpp>
#include <sgpp/combigrid/operation/multidim/fullgrid/AbstractFullGridEvaluator.hpp>
#include <sgpp/combigrid/operation/multidim/fullgrid/AbstractFullGridSummationStrategy.hpp>
#include <sgpp/combigrid/operation/multidim/fullgrid/FullGridLinearSummationStrategy.hpp>
#include <sgpp/combigrid/operation/multidim/fullgrid/FullGridOptimizedPCESummationStrategy.hpp>
#include <sgpp/combigrid/operation/multidim/fullgrid/FullGridPCESummationStrategy.hpp>
#include <sgpp/combigrid/operation/multidim/fullgrid/FullGridQuadraticSummationStrategy.hpp>
#include <sgpp/combigrid/operation/multidim/fullgrid/FullGridTensorVarianceSummationStrategy.hpp>
#include <sgpp/combigrid/operation/multidim/fullgrid/FullGridVarianceSummationStrategy.hpp>
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
typedef GeneralFunction<std::shared_ptr<TreeStorage<double>>, std::shared_ptr<TensorGrid>>
    GridFunction;

template <typename V>
class AbstractFullGridEvaluationStrategy : public AbstractFullGridEvaluator<V> {
 public:
  AbstractFullGridEvaluationStrategy(
      std::shared_ptr<AbstractCombigridStorage> storage,
      std::vector<std::shared_ptr<AbstractLinearEvaluator<V>>> evaluatorPrototypes,
      std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
      FullGridSummationStrategyType summationStrategyType)
      : AbstractFullGridEvaluator<V>(storage, pointHierarchies),
        evaluatorPrototypes(evaluatorPrototypes),
        summationStrategyType(summationStrategyType) {
    switch (summationStrategyType) {
      case FullGridSummationStrategyType::LINEAR:
        summationStrategy = std::make_shared<FullGridLinearSummationStrategy<V>>(
            storage, evaluatorPrototypes, pointHierarchies);
        break;
      case FullGridSummationStrategyType::QUADRATIC:
        summationStrategy = std::make_shared<FullGridQuadraticSummationStrategy<V>>(
            storage, evaluatorPrototypes, pointHierarchies);
        break;
      case FullGridSummationStrategyType::VARIANCE:
        summationStrategy = std::make_shared<FullGridVarianceSummationStrategy<V>>(
            storage, evaluatorPrototypes, pointHierarchies);
        break;
      case FullGridSummationStrategyType::TENSORVARIANCE:
        summationStrategy = std::make_shared<FullGridTensorVarianceSummationStrategy<V>>(
            storage, evaluatorPrototypes, pointHierarchies);
        break;
      case FullGridSummationStrategyType::FULLSUBSPACEDPCE:
        summationStrategy = std::make_shared<FullGridPCESummationStrategy<V>>(
            storage, evaluatorPrototypes, pointHierarchies);
        break;
      case FullGridSummationStrategyType::ONEDSUBSPACEPCE:
        summationStrategy = std::make_shared<FullGridOptimizedPCESummationStrategy<V>>(
            storage, evaluatorPrototypes, pointHierarchies);
        break;
    }
  }

  virtual ~AbstractFullGridEvaluationStrategy() {}

  /**
   * Sets the parameters for the evaluators. Each dimension in which the evaluator does not need a
   * parameter is skipped.
   * So if only the evaluators at dimensions 1 and 3 need a parameter, params.size() should be 2 (or
   * at least 2)
   */
  void setParameters(std::vector<V> const &params) { summationStrategy->setParameters(params); }

  std::vector<std::shared_ptr<AbstractLinearEvaluator<V>>> getEvaluatorPrototypes() {
    return evaluatorPrototypes;
  }

  FullGridSummationStrategyType getSummationStrategyType() { return summationStrategyType; }

 protected:
  std::vector<std::shared_ptr<AbstractLinearEvaluator<V>>> evaluatorPrototypes;
  FullGridSummationStrategyType summationStrategyType;
  std::shared_ptr<AbstractFullGridSummationStrategy<V>> summationStrategy;
};

} /* namespace combigrid */
} /* namespace sgpp */
