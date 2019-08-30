// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/combigrid/algebraic/FloatArrayVector.hpp>
#include <sgpp/combigrid/algebraic/FloatTensorVector.hpp>
#include <sgpp/combigrid/common/MultiIndexIterator.hpp>
#include <sgpp/combigrid/definitions.hpp>
#include <sgpp/combigrid/functions/WeightFunctionsCollection.hpp>
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

// This Summation Strategy calculates the variance on each full grid. The variance is calculated as
//                                  E(u^2) - E(u)^2
// where E(u)^2 is calculated via the scalarProductEvaluatorPrototypes and E(u^2) is calculated via
// a quadrature routine created from the scalarProductEvaluatorPrototypes

template <typename V>
class FullGridVarianceSummationStrategy : public AbstractFullGridSummationStrategy<V> {
 public:
  /**
   * Constructor.
   *
   * @param storage Storage that stores and provides the function values for each grid
   * point.
   * @param scalarProductEvaluatorPrototypes prototype objects for the evaluators that are cloned to
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
                                             pointHierarchies) {
    //    std::cout << "FullGridVarianceSummationStrategy :" << std::endl;
    //    sgpp::combigrid::SingleFunction onedim_weight_function;
    //    for (size_t d = 0; d < scalarProductEvaluatorPrototypes.size(); d++) {
    //      scalarProductEvaluatorPrototypes[d]->getWeightFunction(onedim_weight_function);
    //      std::cout << d << " w(0.5) = " << onedim_weight_function(0.5) << std::endl;
    //  }
  }

  ~FullGridVarianceSummationStrategy() override {}

  /**
   * Calculates the variance by calculating the mean via B spline Quadrature and the mean of f^2 via
   * B spline scalar products
   */
  V eval(MultiIndex const &level) override {
    size_t numDimensions = level.size();

    std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>>
        linearEvaluatorPrototypes;

    EvaluatorConfiguration linearEvalConfig;

    for (size_t d = 0; d < numDimensions; d++) {
      auto evalConfig = this->evaluatorPrototypes[d]->getConfig();
      auto basisType = evalConfig.type;
      if (basisType == CombiEvaluatorTypes::Multi_BSplineScalarProduct) {
        linearEvalConfig.degree = evalConfig.degree;
        linearEvalConfig.type = CombiEvaluatorTypes::Scalar_BSplineQuadrature;
      } else if (basisType == CombiEvaluatorTypes::Multi_PolynomialScalarProduct) {
        linearEvalConfig.type = CombiEvaluatorTypes::Scalar_PolynomialQuadrature;
      } else {
        std::cerr << "FullGridVarianceSummationStrategy: this evaluator is currently "
                     "not supported."
                  << std::endl;
      }
      linearEvaluatorPrototypes.push_back(
          CombiEvaluators::createCombiScalarEvaluator(linearEvalConfig));
    }

    // if a custom weight function shall be used it and its bounds are extracted from the scalar
    // product evaluators

    for (size_t d = 0; d < numDimensions; d++) {
      if (this->evaluatorPrototypes[d]->hasCustomWeightFunction()) {
        sgpp::combigrid::SingleFunction onedim_weight_function;
        double a;
        double b;
        this->evaluatorPrototypes[d]->getWeightFunction(onedim_weight_function);
        this->evaluatorPrototypes[d]->getBounds(a, b);
        linearEvaluatorPrototypes[d]->setWeightFunction(onedim_weight_function);
        linearEvaluatorPrototypes[d]->setBounds(a, b);
      }
    }

    FullGridLinearSummationStrategy<FloatScalarVector> linearStrategy =
        FullGridLinearSummationStrategy<FloatScalarVector>(this->storage, linearEvaluatorPrototypes,
                                                           this->pointHierarchies);

    FloatScalarVector mean = linearStrategy.eval(level);

    // Var = E(u^2) - E(u)^2
    FullGridQuadraticSummationStrategy<V> quadraticStrategyOld =
        FullGridQuadraticSummationStrategy<V>(this->storage, this->evaluatorPrototypes,
                                              this->pointHierarchies);
    V meanSquare = quadraticStrategyOld.eval(level);
    mean.componentwiseMult(mean);
    FloatScalarVector varianceOld = meanSquare[0];
    varianceOld.sub(mean);
    V returnVariance(varianceOld);
    return returnVariance;
  }
};

} /* namespace combigrid */
} /* namespace sgpp */
