// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/operation/Configurations.hpp>

#include <sgpp/combigrid/grid/distribution/ClenshawCurtisDistribution.hpp>
#include <sgpp/combigrid/grid/distribution/LejaPointDistribution.hpp>
#include <sgpp/combigrid/grid/distribution/UniformPointDistribution.hpp>
#include <sgpp/combigrid/grid/growth/ExponentialGrowthStrategy.hpp>
#include <sgpp/combigrid/grid/growth/LinearGrowthStrategy.hpp>
#include <sgpp/combigrid/grid/hierarchy/NestedPointHierarchy.hpp>
#include <sgpp/combigrid/grid/hierarchy/NonNestedPointHierarchy.hpp>
#include <sgpp/combigrid/grid/ordering/ExponentialLevelorderPointOrdering.hpp>
#include <sgpp/combigrid/grid/ordering/IdentityPointOrdering.hpp>

#include <sgpp/combigrid/operation/onedim/ArrayEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/LinearInterpolationEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/PolynomialInterpolationEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/QuadratureEvaluator.hpp>

namespace sgpp {
namespace combigrid {

std::shared_ptr<AbstractPointHierarchy> CombiHierarchies::linearLeja(size_t growthFactor) {
  return std::make_shared<NestedPointHierarchy>(
      std::make_shared<LejaPointDistribution>(),
      std::make_shared<IdentityPointOrdering>(std::make_shared<LinearGrowthStrategy>(growthFactor),
                                              false));
}

std::shared_ptr<AbstractPointHierarchy> CombiHierarchies::expLeja() {
  return std::make_shared<NestedPointHierarchy>(
      std::make_shared<LejaPointDistribution>(),
      std::make_shared<IdentityPointOrdering>(std::make_shared<ExponentialGrowthStrategy>(),
                                              false));
}

std::shared_ptr<AbstractPointHierarchy> CombiHierarchies::expUniform() {
  return std::make_shared<NestedPointHierarchy>(
      std::make_shared<UniformPointDistribution>(),
      std::make_shared<ExponentialLevelorderPointOrdering>());
}

std::shared_ptr<AbstractPointHierarchy> CombiHierarchies::expClenshawCurtis() {
  return std::make_shared<NestedPointHierarchy>(
      std::make_shared<ClenshawCurtisDistribution>(),
      std::make_shared<ExponentialLevelorderPointOrdering>());
}

std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>
CombiEvaluators::polynomialInterpolation() {
  return std::make_shared<PolynomialInterpolationEvaluator>();
}

std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>> CombiEvaluators::linearInterpolation() {
  return std::make_shared<LinearInterpolationEvaluator>();
}

std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>> CombiEvaluators::quadrature() {
  return std::make_shared<QuadratureEvaluator>();
}

std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>>
CombiEvaluators::multiPolynomialInterpolation() {
  return std::make_shared<ArrayEvaluator<PolynomialInterpolationEvaluator>>(true);
}

std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>>
CombiEvaluators::multiLinearInterpolation() {
  return std::make_shared<ArrayEvaluator<LinearInterpolationEvaluator>>(true);
}

std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>> CombiEvaluators::multiQuadrature() {
  return std::make_shared<ArrayEvaluator<QuadratureEvaluator>>(false);
}

} /* namespace combigrid */
} /* namespace sgpp */
