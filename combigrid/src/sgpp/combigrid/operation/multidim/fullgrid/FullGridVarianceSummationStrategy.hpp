// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/combigrid/algebraic/FloatArrayVector.hpp>
#include <sgpp/combigrid/algebraic/FloatTensorVector.hpp>
#include <sgpp/combigrid/common/MultiIndexIterator.hpp>
#include <sgpp/combigrid/definitions.hpp>
#include <sgpp/combigrid/grid/hierarchy/AbstractPointHierarchy.hpp>
#include <sgpp/combigrid/operation/multidim/fullgrid/AbstractFullGridSummationStrategy.hpp>
#include <sgpp/combigrid/operation/multidim/fullgrid/FullGridLinearSummationStrategy.hpp>
#include <sgpp/combigrid/operation/multidim/fullgrid/FullGridQuadraticSummationStrategy.hpp>
#include <sgpp/combigrid/storage/AbstractCombigridStorage.hpp>
#include <sgpp/combigrid/threading/PtrGuard.hpp>
#include <sgpp/combigrid/threading/ThreadPool.hpp>

#include <iostream>
#include <vector>

namespace sgpp {
namespace combigrid {

template <typename V>
class FullGridVarianceSummationStrategy : public AbstractFullGridSummationStrategy<V> {
 protected:
 private:
  std::shared_ptr<AbstractCombigridStorage> storage;
  std::vector<std::shared_ptr<AbstractLinearEvaluator<V>>> linearEvaluatorPrototypes;
  std::vector<std::shared_ptr<AbstractLinearEvaluator<V>>> quadraticEvaluatorPrototypes;
  std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies;
  std::vector<std::vector<V>> linearBasisValues;
  std::vector<std::vector<V>> quadraticBasisValues;

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
  FullGridVarianceSummationStrategy(
      std::shared_ptr<AbstractCombigridStorage> storage,
      std::vector<std::shared_ptr<AbstractLinearEvaluator<V>>> linearEvaluatorPrototypes,
      std::vector<std::shared_ptr<AbstractLinearEvaluator<V>>> quadraticEvaluatorPrototypes,
      std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies)
      : storage(storage),
        linearEvaluatorPrototypes(linearEvaluatorPrototypes),
        quadraticEvaluatorPrototypes(quadraticEvaluatorPrototypes),
        pointHierarchies(pointHierarchies),
        linearBasisValues(linearEvaluatorPrototypes.size()),
        quadraticBasisValues(quadraticEvaluatorPrototypes.size()) {}

  ~FullGridVarianceSummationStrategy() {}

  /**
   * SOME DESCRIPTION
   */
  V eval(MultiIndex const &level) override {
    auto linearStrategy = FullGridLinearSummationStrategy<FloatScalarVector>(
        storage, linearEvaluatorPrototypes, pointHierarchies);
    auto quadraticStrategy = FullGridLinearSummationStrategy<FloatArrayVector>(
        this->storage, quadraticEvaluatorPrototypes, this->pointHierarchies);

    FloatScalarVector mean = linearStrategy.eval(level);
    FloatArrayVector meanSquare = quadraticStrategy.eval(level);

    // Var = E(x^2) - E(x)^2
    mean.componentwiseMult(mean);
    meanSquare[0].sub(mean);
    FloatScalarVector variance = meanSquare[0].value();

    V returnVariance(variance);
    return returnVariance;
  }
};

} /* namespace combigrid */
} /* namespace sgpp */
