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
#include <sgpp/combigrid/operation/Configurations.hpp>
#include <sgpp/combigrid/operation/multidim/fullgrid/AbstractFullGridSummationStrategy.hpp>
#include <sgpp/combigrid/operation/multidim/fullgrid/FullGridLinearSummationStrategy.hpp>
#include <sgpp/combigrid/operation/multidim/fullgrid/FullGridQuadraticSummationStrategy.hpp>
#include <sgpp/combigrid/operation/onedim/AbstractLinearEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/BSplineQuadratureEvaluator.hpp>
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
  //  std::shared_ptr<AbstractCombigridStorage> storage;
  //  std::vector<std::shared_ptr<AbstractLinearEvaluator<V>>> scalarProductEvaluatorPrototypes;
  //  std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies;

 public:
  /**
   * Constructor.
   *
   * @param storage Storage that stores and provides the function values for each grid
   * point.
   * @param evaluatorPrototypes prototype objects for the evaluators that are cloned to
   * get an
   * evaluator for each dimension and each level.
   * @param pointHierarchies PointHierarchy objects for each dimension providing the
   * points for each
   * level and information about their ordering.
   */
  FullGridVarianceSummationStrategy(
      std::shared_ptr<AbstractCombigridStorage> storage,
      std::vector<std::shared_ptr<AbstractLinearEvaluator<V>>> scalarProductEvaluatorPrototypes,
      std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies)
      : AbstractFullGridSummationStrategy<V>(storage, scalarProductEvaluatorPrototypes,
                                             pointHierarchies) {}

  ~FullGridVarianceSummationStrategy() {}

  /**
   * SOME DESCRIPTION
   * People want to know everything about this method including the secret Bspline
   * techniques!
   *
   * Currently only V=FloatArrayVector is supported
   */
  V eval(MultiIndex const &level) override {
    auto evalConfig = this->evaluatorPrototypes[0]->getConfig();
    size_t numDimensions = level.size();

    bool bsplineBasis = true;
    for (size_t d = 0; d < numDimensions; d++) {
      if (evalConfig.type != CombiEvaluatorTypes::Multi_BSplineScalarProduct) {
        bsplineBasis = false;
        std::cerr << "FullGridVarianceSummationStrategy: this evaluator is currently "
                     "not supported."
                  << std::endl;
      }
    }

    if (bsplineBasis) {
      FullGridQuadraticSummationStrategy<V> quadraticStrategy =
          FullGridQuadraticSummationStrategy<V>(this->storage, this->evaluatorPrototypes,
                                                this->pointHierarchies);

      size_t degree = evalConfig.degree;
      EvaluatorConfiguration linearEvalConfig(CombiEvaluatorTypes::Scalar_BSplineQuadrature,
                                              degree);
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>>
          linearEvaluatorPrototypes(this->evaluatorPrototypes.size(),
                                    CombiEvaluators::createCombiScalarEvaluator(linearEvalConfig));

      FullGridLinearSummationStrategy<FloatScalarVector> linearStrategy =
          FullGridLinearSummationStrategy<FloatScalarVector>(
              this->storage, linearEvaluatorPrototypes, this->pointHierarchies);

      V meanSquare = quadraticStrategy.eval(level);
      FloatScalarVector mean = linearStrategy.eval(level);
      std::cout.precision(10);
      //      std::cout << "mean " << mean.value() << " meanSquare " << meanSquare[0].value();

      // Var = E(x^2) - E(x)^2
      mean.componentwiseMult(mean);
      FloatScalarVector variance = meanSquare[0];
      variance.sub(mean);

      //      std::cout << " variance " << variance.value() << std::endl;

      V returnVariance(variance);
      return returnVariance;
    }
    V returnFailure = V::zero();
    std::cerr << "FullGridVarianceSummationStrategy: this evaluator is currently "
                 "not supported. Returning zero."
              << std::endl;
    return returnFailure;
  }
};

} /* namespace combigrid */
} /* namespace sgpp */
