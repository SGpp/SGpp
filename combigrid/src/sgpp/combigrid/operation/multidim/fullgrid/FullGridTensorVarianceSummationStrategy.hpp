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

#include <sgpp/base/exception/algorithm_exception.hpp>

#include <iostream>
#include <vector>

namespace sgpp {
namespace combigrid {

template <typename V>
class FullGridTensorVarianceSummationStrategy : public AbstractFullGridSummationStrategy<V> {
 public:
  /**
   * Constructor.
   *
   * @param storage Storage that stores and provides the function values for each grid
   * point.
   * @param scalarProductEvaluatorPrototypes prototype objects for the evaluators that are cloned to
   * get an evaluator for each dimension and each level.
   * @param pointHierarchies PointHierarchy objects for each dimension providing the
   * points for each
   * level and information about their ordering.
   */
  FullGridTensorVarianceSummationStrategy(
      std::shared_ptr<AbstractCombigridStorage> storage,
      std::vector<std::shared_ptr<AbstractLinearEvaluator<V>>> scalarProductEvaluatorPrototypes,
      std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies)
      : AbstractFullGridSummationStrategy<V>(storage, scalarProductEvaluatorPrototypes,
                                             pointHierarchies) {}

  ~FullGridTensorVarianceSummationStrategy() override {}

  V eval(MultiIndex const &level) override {
    size_t numDimensions = level.size();

    std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatTensorVector>>>
        linearEvaluatorPrototypes;

    EvaluatorConfiguration linearEvalConfig;

    for (size_t d = 0; d < numDimensions; d++) {
      auto evalConfig = this->evaluatorPrototypes[d]->getConfig();
      auto basisType = evalConfig.type;
      if (basisType == CombiEvaluatorTypes::Tensor_PolynomialInterpolation) {
        linearEvalConfig.type = CombiEvaluatorTypes::Tensor_PolynomialInterpolation;
        linearEvalConfig.functionBasis = evalConfig.functionBasis;
      } else {
        throw sgpp::base::algorithm_exception(
            "FullGridTensorNormSummationStrategy: this evaluator is currently "
            "not supported.");
      }
      linearEvaluatorPrototypes.push_back(
          CombiEvaluators::createCombiTensorEvaluator(linearEvalConfig));
    }

    FullGridLinearSummationStrategy<FloatTensorVector> linearStrategy =
        FullGridLinearSummationStrategy<FloatTensorVector>(this->storage, linearEvaluatorPrototypes,
                                                           this->pointHierarchies);

    // -------------------------------------------------------------------
    // compute the norm of the results and skip the constant part. This
    // is equivalent to the variance if the function basis is an orthogonal one
    // with respect to the marginal densities
    FloatTensorVector result = linearStrategy.eval(level);
    double var = 0.0;
    auto it = result.getValues()->getStoredDataIterator();
    if (it->isValid()) {
      it->moveToNext();  // ignore first entry (belonging to mean)

      for (; it->isValid(); it->moveToNext()) {
        double coeff = it->value().value();
        var += coeff * coeff;
      }
    }
    // -----------------------------------------------------------
    FloatScalarVector var_vec(var);
    V returnVar(var_vec);
    return returnVar;
  }
};

} /* namespace combigrid */
} /* namespace sgpp */
