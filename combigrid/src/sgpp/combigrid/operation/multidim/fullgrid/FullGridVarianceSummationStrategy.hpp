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
                                             pointHierarchies) {}

  ~FullGridVarianceSummationStrategy() {}

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
    double width = 1.0;
    for (size_t d = 0; d < numDimensions; d++) {
      if (this->evaluatorPrototypes[d]->hasCustomWeightFunction()) {
        sgpp::combigrid::SingleFunction weight_function;
        double a;
        double b;
        this->evaluatorPrototypes[d]->getWeightFunction(weight_function);
        this->evaluatorPrototypes[d]->getBounds(a, b);
        linearEvaluatorPrototypes[d]->setWeightFunction(weight_function);
        linearEvaluatorPrototypes[d]->setBounds(a, b);
        width *= b - a;
      }
    }

    FullGridLinearSummationStrategy<FloatScalarVector> linearStrategy =
        FullGridLinearSummationStrategy<FloatScalarVector>(this->storage, linearEvaluatorPrototypes,
                                                           this->pointHierarchies);

    FullGridQuadraticSummationStrategy<V> quadraticStrategy = FullGridQuadraticSummationStrategy<V>(
        this->storage, this->evaluatorPrototypes, this->pointHierarchies);

    // ToDo (rehmemk) first quadrature => numAdditionalPoints for scalar products?

    // Var = E(x^2) - E(x)^2
    FloatScalarVector mean = linearStrategy.eval(level);
    V meanSquare = quadraticStrategy.eval(level);
    mean.scalarMult(width);

    mean.componentwiseMult(mean);
    FloatScalarVector variance = meanSquare[0];

    //    double exactMeanSquare = 0.237373737373740;
    //    std::cout << "(" << level[0] << " " << level[1] << ")  meanSquare : " <<
    //    meanSquare[0].value()
    //              << " error: " << std::fabs(exactMeanSquare - meanSquare[0].value()) <<
    //              std::endl;

    variance.scalarMult(width);
    variance.sub(mean);

    //    double exactVariance = 0.126262626262629;
    //    std::cout << "(" << level[0] << " " << level[1] << ")  vaiance : " << variance[0]
    //                  << " error: " << std::fabs(exactVaraince - variance[0].value()) <<
    //                  std::endl;

    V returnVariance(variance);
    return returnVariance;
  }
};

} /* namespace combigrid */
} /* namespace sgpp */
