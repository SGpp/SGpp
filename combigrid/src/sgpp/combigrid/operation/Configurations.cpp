// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/grid/distribution/ChebyshevDistribution.hpp>
#include <sgpp/combigrid/grid/distribution/ClenshawCurtisDistribution.hpp>
#include <sgpp/combigrid/grid/distribution/L2LejaPointDistribution.hpp>
#include <sgpp/combigrid/grid/distribution/LejaPointDistribution.hpp>
#include <sgpp/combigrid/grid/distribution/UniformBoundaryPointDistribution.hpp>
#include <sgpp/combigrid/grid/distribution/UniformPointDistribution.hpp>
#include <sgpp/combigrid/grid/growth/ExponentialGrowthStrategy.hpp>
#include <sgpp/combigrid/grid/growth/LinearGrowthStrategy.hpp>
#include <sgpp/combigrid/grid/hierarchy/NestedPointHierarchy.hpp>
#include <sgpp/combigrid/grid/hierarchy/NonNestedPointHierarchy.hpp>
#include <sgpp/combigrid/grid/ordering/ExponentialChebyshevPointOrdering.hpp>
#include <sgpp/combigrid/grid/ordering/ExponentialLevelorderPointOrdering.hpp>
#include <sgpp/combigrid/grid/ordering/ExponentialNoBoundaryPointOrdering.hpp>
#include <sgpp/combigrid/grid/ordering/IdentityPointOrdering.hpp>
#include <sgpp/combigrid/operation/Configurations.hpp>

#include <sgpp/base/exception/factory_exception.hpp>
#include <sgpp/combigrid/operation/OperationConfiguration.hpp>
#include <sgpp/combigrid/operation/onedim/ArrayEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/BSplineInterpolationCoefficientEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/BSplineInterpolationEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/BSplineQuadratureEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/BSplineScalarProductEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/CubicSplineInterpolationEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/InterpolationCoefficientEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/LinearInterpolationEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/PolynomialInterpolationEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/PolynomialQuadratureEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/PolynomialScalarProductEvaluator.hpp>

namespace sgpp {
namespace combigrid {

std::shared_ptr<AbstractPointHierarchy> CombiHierarchies::linearLeja(size_t growthFactor) {
  return std::make_shared<NestedPointHierarchy>(
      std::make_shared<LejaPointDistribution>(),
      std::make_shared<IdentityPointOrdering>(std::make_shared<LinearGrowthStrategy>(growthFactor),
                                              false));
}

std::shared_ptr<AbstractPointHierarchy> CombiHierarchies::linearLeja(SingleFunction weightFunction,
                                                                     size_t growthFactor) {
  return std::make_shared<NestedPointHierarchy>(
      std::make_shared<LejaPointDistribution>(weightFunction),
      std::make_shared<IdentityPointOrdering>(std::make_shared<LinearGrowthStrategy>(growthFactor),
                                              false));
}

std::shared_ptr<AbstractPointHierarchy> CombiHierarchies::linearL2Leja(size_t growthFactor) {
  return std::make_shared<NestedPointHierarchy>(
      std::make_shared<L2LejaPointDistribution>(),
      std::make_shared<IdentityPointOrdering>(std::make_shared<LinearGrowthStrategy>(growthFactor),
                                              false));
}

std::shared_ptr<AbstractPointHierarchy> CombiHierarchies::linearL2Leja(
    SingleFunction weightFunction, size_t growthFactor, size_t numAdditionalPoints) {
  return std::make_shared<NestedPointHierarchy>(
      std::make_shared<L2LejaPointDistribution>(weightFunction, numAdditionalPoints),
      std::make_shared<IdentityPointOrdering>(std::make_shared<LinearGrowthStrategy>(growthFactor),
                                              false));
}

std::shared_ptr<AbstractPointHierarchy> CombiHierarchies::expLeja() {
  return std::make_shared<NestedPointHierarchy>(
      std::make_shared<LejaPointDistribution>(),
      std::make_shared<IdentityPointOrdering>(std::make_shared<ExponentialGrowthStrategy>(),
                                              false));
}

std::shared_ptr<AbstractPointHierarchy> CombiHierarchies::expLeja(SingleFunction weightFunction) {
  return std::make_shared<NestedPointHierarchy>(
      std::make_shared<LejaPointDistribution>(weightFunction),
      std::make_shared<IdentityPointOrdering>(std::make_shared<ExponentialGrowthStrategy>(),
                                              false));
}

std::shared_ptr<AbstractPointHierarchy> CombiHierarchies::expL2Leja() {
  return std::make_shared<NestedPointHierarchy>(
      std::make_shared<L2LejaPointDistribution>(),
      std::make_shared<IdentityPointOrdering>(std::make_shared<ExponentialGrowthStrategy>(),
                                              false));
}

std::shared_ptr<AbstractPointHierarchy> CombiHierarchies::expL2Leja(SingleFunction weightFunction,
                                                                    size_t numAdditionalPoints) {
  return std::make_shared<NestedPointHierarchy>(
      std::make_shared<L2LejaPointDistribution>(weightFunction, numAdditionalPoints),
      std::make_shared<IdentityPointOrdering>(std::make_shared<ExponentialGrowthStrategy>(),
                                              false));
}

std::shared_ptr<AbstractPointHierarchy> CombiHierarchies::expUniform() {
  return std::make_shared<NestedPointHierarchy>(
      std::make_shared<UniformPointDistribution>(),
      std::make_shared<ExponentialNoBoundaryPointOrdering>());
}

std::shared_ptr<AbstractPointHierarchy> CombiHierarchies::expClenshawCurtis() {
  return std::make_shared<NestedPointHierarchy>(
      std::make_shared<ClenshawCurtisDistribution>(),
      std::make_shared<ExponentialLevelorderPointOrdering>());
}

std::shared_ptr<AbstractPointHierarchy> CombiHierarchies::expChebyshev() {
  return std::make_shared<NestedPointHierarchy>(
      std::make_shared<ChebyshevDistribution>(),
      std::make_shared<ExponentialChebyshevPointOrdering>());
}

std::shared_ptr<AbstractPointHierarchy> CombiHierarchies::expUniformBoundary() {
  return std::make_shared<NestedPointHierarchy>(
      std::make_shared<UniformBoundaryPointDistribution>(),
      std::make_shared<ExponentialLevelorderPointOrdering>());
}

std::shared_ptr<AbstractPointHierarchy> CombiHierarchies::linearUniform(size_t growthFactor) {
  return std::make_shared<NonNestedPointHierarchy>(
      std::make_shared<UniformPointDistribution>(),
      std::make_shared<IdentityPointOrdering>(std::make_shared<LinearGrowthStrategy>(growthFactor),
                                              true));
}

std::shared_ptr<AbstractPointHierarchy> CombiHierarchies::linearClenshawCurtis(
    size_t growthFactor) {
  return std::make_shared<NonNestedPointHierarchy>(
      std::make_shared<ClenshawCurtisDistribution>(),
      std::make_shared<IdentityPointOrdering>(std::make_shared<LinearGrowthStrategy>(growthFactor),
                                              true));
}

std::shared_ptr<AbstractPointHierarchy> CombiHierarchies::linearChebyshev(size_t growthFactor) {
  return std::make_shared<NonNestedPointHierarchy>(
      std::make_shared<ChebyshevDistribution>(),
      std::make_shared<IdentityPointOrdering>(std::make_shared<LinearGrowthStrategy>(growthFactor),
                                              true));
}

std::shared_ptr<AbstractPointHierarchy> CombiHierarchies::linearUniformBoundary(
    size_t growthFactor) {
  return std::make_shared<NonNestedPointHierarchy>(
      std::make_shared<UniformBoundaryPointDistribution>(),
      std::make_shared<IdentityPointOrdering>(std::make_shared<LinearGrowthStrategy>(growthFactor),
                                              true));
}

std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>
CombiEvaluators::createCombiScalarEvaluator(EvaluatorConfiguration evalConfig) {
  if (evalConfig.type == CombiEvaluatorTypes::Scalar_PolynomialInterpolation) {
    return std::make_shared<PolynomialInterpolationEvaluator>();
  }
  if (evalConfig.type == CombiEvaluatorTypes::Scalar_LinearInterpolation) {
    return std::make_shared<LinearInterpolationEvaluator>();
  }
  if (evalConfig.type == CombiEvaluatorTypes::Scalar_CubicSplineInterpolation) {
    return std::make_shared<CubicSplineInterpolationEvaluator>();
  }
  if (evalConfig.type == CombiEvaluatorTypes::Scalar_BSplineInterpolation) {
    return std::make_shared<BSplineInterpolationEvaluator>(evalConfig.degree);
  }
  if (evalConfig.type == CombiEvaluatorTypes::Scalar_PolynomialQuadrature) {
    return std::make_shared<PolynomialQuadratureEvaluator>();
  }
  if (evalConfig.type == CombiEvaluatorTypes::Scalar_BSplineQuadrature) {
    return std::make_shared<BSplineQuadratureEvaluator>(evalConfig.degree);
  } else {
    throw sgpp::base::factory_exception(
        "CombiEvaluators::createCombiScalarEvaluator: type is not supported");
  }
}

std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>
CombiEvaluators::polynomialInterpolation() {
  return std::make_shared<PolynomialInterpolationEvaluator>();
}

std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>> CombiEvaluators::linearInterpolation() {
  return std::make_shared<LinearInterpolationEvaluator>();
}

std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>
CombiEvaluators::cubicSplineInterpolation() {
  return std::make_shared<CubicSplineInterpolationEvaluator>();
}

std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>> CombiEvaluators::BSplineInterpolation(
    size_t degree) {
  return std::make_shared<BSplineInterpolationEvaluator>(degree);
}

std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>> CombiEvaluators::quadrature() {
  return std::make_shared<PolynomialQuadratureEvaluator>();
}

std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>> CombiEvaluators::BSplineQuadrature(
    size_t degree) {
  return std::make_shared<BSplineQuadratureEvaluator>(degree);
}

std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>> CombiEvaluators::BSplineQuadrature(
    size_t degree, sgpp::combigrid::SingleFunction weight_function, size_t numAdditionalPoints,
    bool normalizeWeights) {
  return std::make_shared<BSplineQuadratureEvaluator>(degree, weight_function, numAdditionalPoints,
                                                      normalizeWeights);
}

std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>> CombiEvaluators::BSplineQuadrature(
    size_t degree, sgpp::combigrid::SingleFunction weight_function, size_t numAdditionalPoints,
    double a, double b, bool normalizeWeights) {
  return std::make_shared<BSplineQuadratureEvaluator>(degree, weight_function, numAdditionalPoints,
                                                      a, b, normalizeWeights);
}

std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>>
CombiEvaluators::createCombiMultiEvaluator(EvaluatorConfiguration evalConfig) {
  if (evalConfig.type == CombiEvaluatorTypes::Multi_PolynomialInterpolation) {
    return std::make_shared<ArrayEvaluator<PolynomialInterpolationEvaluator>>(true);
  } else if (evalConfig.type == CombiEvaluatorTypes::Multi_PolynomialQuadrature) {
    return std::make_shared<ArrayEvaluator<PolynomialQuadratureEvaluator>>(false);
  } else if (evalConfig.type == CombiEvaluatorTypes::Multi_BSplineInterpolation) {
    return std::make_shared<ArrayEvaluator<BSplineInterpolationEvaluator>>(
        true, BSplineInterpolationEvaluator(evalConfig.degree));
  } else if (evalConfig.type == CombiEvaluatorTypes::Multi_BSplineQuadrature) {
    return std::make_shared<ArrayEvaluator<BSplineQuadratureEvaluator>>(
        false, BSplineQuadratureEvaluator(evalConfig.degree));
  } else if (evalConfig.type == CombiEvaluatorTypes::Multi_BSplineScalarProduct) {
    std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>> evaluator =
        std::make_shared<BSplineScalarProductEvaluator>(evalConfig.degree);
    evaluator->setConfig(evalConfig);
    return evaluator;
  } else if (evalConfig.type == CombiEvaluatorTypes::Multi_PolynomialScalarProduct) {
    std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>> evaluator =
        std::make_shared<PolynomialScalarProductEvaluator>();
    evaluator->setConfig(evalConfig);
    return evaluator;
  } else if (evalConfig.type == CombiEvaluatorTypes::Multi_CubicSplineInterpolation) {
    return std::make_shared<ArrayEvaluator<CubicSplineInterpolationEvaluator>>(true);
  } else if (evalConfig.type == CombiEvaluatorTypes::Multi_LinearInterpolation) {
    return std::make_shared<ArrayEvaluator<LinearInterpolationEvaluator>>(true);
  } else {
    throw sgpp::base::factory_exception(
        "CombiEvaluators::createCombiMultiEvaluator: type is not supported");
  }
}

std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>> CombiEvaluators::BSplineScalarProduct(
    size_t degree) {
  return std::make_shared<BSplineScalarProductEvaluator>(degree);
}

std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>>
CombiEvaluators::polynomialScalarProduct() {
  return std::make_shared<PolynomialScalarProductEvaluator>();
}

std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>>
CombiEvaluators::multiPolynomialInterpolation() {
  return std::make_shared<ArrayEvaluator<PolynomialInterpolationEvaluator>>(true);
}

std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>>
CombiEvaluators::multiLinearInterpolation() {
  return std::make_shared<ArrayEvaluator<LinearInterpolationEvaluator>>(true);
}

std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>>
CombiEvaluators::multiCubicSplineInterpolation() {
  return std::make_shared<ArrayEvaluator<CubicSplineInterpolationEvaluator>>(true);
}

std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>>
CombiEvaluators::multiBSplineInterpolation(size_t degree) {
  return std::make_shared<ArrayEvaluator<BSplineInterpolationEvaluator>>(
      true, BSplineInterpolationEvaluator(degree));
}

std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>> CombiEvaluators::multiBSplineQuadrature(
    size_t degree) {
  return std::make_shared<ArrayEvaluator<BSplineQuadratureEvaluator>>(
      false, BSplineQuadratureEvaluator(degree));
}

std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>> CombiEvaluators::multiQuadrature() {
  return std::make_shared<ArrayEvaluator<PolynomialQuadratureEvaluator>>(false);
}

std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>> CombiEvaluators::multiQuadrature(
    SingleFunction func, bool normalizeWeights) {
  return std::make_shared<ArrayEvaluator<PolynomialQuadratureEvaluator>>(
      false, PolynomialQuadratureEvaluator(func, normalizeWeights));
}

// tensor evaluators
std::shared_ptr<AbstractLinearEvaluator<FloatTensorVector>>
CombiEvaluators::createCombiTensorEvaluator(EvaluatorConfiguration evalConfig) {
  if (evalConfig.type == CombiEvaluatorTypes::Tensor_PolynomialInterpolation) {
    return std::make_shared<InterpolationCoefficientEvaluator>(evalConfig.functionBasis);
  } else {
    throw sgpp::base::factory_exception(
        "CombiEvaluators::createCombiTensorEvaluator: type is not supported");
  }
}

std::shared_ptr<AbstractLinearEvaluator<FloatTensorVector>> CombiEvaluators::tensorInterpolation(
    std::shared_ptr<AbstractInfiniteFunctionBasis1D> functionBasis) {
  return std::make_shared<InterpolationCoefficientEvaluator>(functionBasis);
}

std::shared_ptr<AbstractLinearEvaluator<FloatTensorVector>>
CombiEvaluators::tensorBSplineInterpolation(size_t degree) {
  return std::make_shared<BSplineInterpolationCoefficientEvaluator>(degree);
}

} /* namespace combigrid */
} /* namespace sgpp */
