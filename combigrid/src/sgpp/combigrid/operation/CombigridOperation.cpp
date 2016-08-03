// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/grid/distribution/ClenshawCurtisDistribution.hpp>
#include <sgpp/combigrid/grid/growth/ExponentialGrowthStrategy.hpp>
#include <sgpp/combigrid/grid/hierarchy/NestedPointHierarchy.hpp>
#include <sgpp/combigrid/grid/ordering/ExponentialLevelorderPointOrdering.hpp>
#include <sgpp/combigrid/operation/CombigridOperation.hpp>
#include <sgpp/combigrid/operation/Configurations.hpp>
#include <sgpp/combigrid/operation/multidim/CombigridEvaluator.hpp>
#include <sgpp/combigrid/operation/multidim/FullGridTensorEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/BarycentricInterpolationEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/LinearInterpolationEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/QuadratureEvaluator.hpp>
#include <sgpp/combigrid/storage/tree/CombigridTreeStorage.hpp>

#include <vector>

namespace sgpp {
namespace combigrid {

class CombigridOperationImpl {
 public:
  CombigridOperationImpl(
      std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>> evaluatorPrototypes,
      std::shared_ptr<AbstractCombigridStorage> storage)
      : storage(storage),
        fullGridEval(new FullGridTensorEvaluator<FloatScalarVector>(storage, evaluatorPrototypes,
                                                                    pointHierarchies)),
        combiEval(
            new CombigridEvaluator<FloatScalarVector>(pointHierarchies.size(), fullGridEval)) {}

  std::shared_ptr<AbstractCombigridStorage> storage;
  std::shared_ptr<FullGridTensorEvaluator<FloatScalarVector>> fullGridEval;
  std::shared_ptr<CombigridEvaluator<FloatScalarVector>> combiEval;
};

CombigridOperation::CombigridOperation(
    std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
    std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>> evaluatorPrototypes,
    MultiFunction func)
    : impl(new CombigridOperationImpl(pointHierarchies, evaluatorPrototypes,
                                      std::shared_ptr<AbstractCombigridStorage>(
                                          new CombigridTreeStorage(pointHierarchies, func)))) {}

CombigridOperation::CombigridOperation(
    std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
    std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>> evaluatorPrototypes,
    std::shared_ptr<AbstractCombigridStorage> storage)
    : impl(new CombigridOperationImpl(pointHierarchies, evaluatorPrototypes, storage)) {}

void CombigridOperation::setParameters(const base::DataVector& param) {
  std::vector<FloatScalarVector> scalars(param.getSize());
  for (size_t i = 0; i < param.getSize(); ++i) {
    scalars[i].value() = param[i];
  }

  impl->fullGridEval->setParameters(scalars);
  impl->combiEval->clear();
}

double CombigridOperation::getResult() { return impl->combiEval->getValue().value(); }

double CombigridOperation::evaluate(size_t q, base::DataVector const& param) {
  setParameters(param);
  impl->combiEval->addRegularLevels(q);

  return getResult();
}

std::shared_ptr<CombigridOperation>
CombigridOperation::createExpClenshawCurtisPolynomialInterpolation(size_t numDimensions,
                                                                   MultiFunction func) {
  return std::make_shared<CombigridOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(numDimensions,
                                                           CombiHierarchies::expClenshawCurtis()),
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>>(
          numDimensions, CombiEvaluators::polynomialInterpolation()),
      func);
}

std::shared_ptr<CombigridOperation> CombigridOperation::createExpLejaPolynomialInterpolation(
    size_t numDimensions, MultiFunction func) {
  return std::make_shared<CombigridOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(numDimensions,
                                                           CombiHierarchies::expLeja()),
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>>(
          numDimensions, CombiEvaluators::polynomialInterpolation()),
      func);
}

std::shared_ptr<CombigridOperation> CombigridOperation::createExpUniformPolynomialInterpolation(
    size_t numDimensions, MultiFunction func) {
  return std::make_shared<CombigridOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(numDimensions,
                                                           CombiHierarchies::expUniform()),
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>>(
          numDimensions, CombiEvaluators::polynomialInterpolation()),
      func);
}

std::shared_ptr<CombigridOperation>
CombigridOperation::createLinearClenshawCurtisPolynomialInterpolation(size_t numDimensions,
                                                                      MultiFunction func) {
  return std::make_shared<CombigridOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(
          numDimensions, std::make_shared<NonNestedPointHierarchy>(
                             std::make_shared<ClenshawCurtisDistribution>(),
                             std::make_shared<IdentityPointOrdering>(
                                 std::make_shared<LinearGrowthStrategy>(2), true))),
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>>(
          numDimensions, CombiEvaluators::polynomialInterpolation()),
      func);
}

std::shared_ptr<CombigridOperation> CombigridOperation::createLinearLejaPolynomialInterpolation(
    size_t numDimensions, MultiFunction func, size_t growthFactor) {
  return std::make_shared<CombigridOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(
          numDimensions, CombiHierarchies::linearLeja(growthFactor)),
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>>(
          numDimensions, CombiEvaluators::polynomialInterpolation()),
      func);
}

std::shared_ptr<CombigridOperation> CombigridOperation::createLinearUniformPolynomialInterpolation(
    size_t numDimensions, MultiFunction func) {
  return std::make_shared<CombigridOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(
          numDimensions, std::make_shared<NonNestedPointHierarchy>(
                             std::make_shared<UniformPointDistribution>(),
                             std::make_shared<IdentityPointOrdering>(
                                 std::make_shared<LinearGrowthStrategy>(2), true))),
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>>(
          numDimensions, CombiEvaluators::polynomialInterpolation()),
      func);
}

std::shared_ptr<CombigridOperation> CombigridOperation::createLinearLejaQuadrature(
    size_t numDimensions, MultiFunction func, size_t growthFactor) {
  return std::make_shared<CombigridOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(
          numDimensions, CombiHierarchies::linearLeja(growthFactor)),
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>>(
          numDimensions, CombiEvaluators::quadrature()),
      func);
}

std::shared_ptr<CombigridOperation> CombigridOperation::createExpUniformLinearInterpolation(
    size_t numDimensions, MultiFunction func) {
  return std::make_shared<CombigridOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(numDimensions,
                                                           CombiHierarchies::expUniform()),
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>>(
          numDimensions, CombiEvaluators::linearInterpolation()),
      func);
}
} /* namespace combigrid */
} /* namespace sgpp*/
