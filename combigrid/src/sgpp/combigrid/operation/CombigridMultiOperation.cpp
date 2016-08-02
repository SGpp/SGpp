// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/operation/CombigridMultiOperation.hpp>
#include <sgpp/combigrid/operation/multidim/CombigridEvaluator.hpp>
#include <sgpp/combigrid/operation/multidim/FullGridTensorEvaluator.hpp>
#include <sgpp/combigrid/storage/tree/CombigridTreeStorage.hpp>
#include <sgpp/combigrid/grid/hierarchy/NestedPointHierarchy.hpp>
#include <sgpp/combigrid/grid/hierarchy/NonNestedPointHierarchy.hpp>
#include <sgpp/combigrid/grid/ordering/ExponentialLevelorderPointOrdering.hpp>
#include <sgpp/combigrid/grid/ordering/IdentityPointOrdering.hpp>
#include <sgpp/combigrid/grid/growth/LinearGrowthStrategy.hpp>
#include <sgpp/combigrid/grid/growth/ExponentialGrowthStrategy.hpp>
#include <sgpp/combigrid/grid/distribution/ClenshawCurtisDistribution.hpp>
#include <sgpp/combigrid/grid/distribution/UniformPointDistribution.hpp>
#include <sgpp/combigrid/grid/distribution/LejaPointDistribution.hpp>
#include <sgpp/combigrid/operation/onedim/BarycentricInterpolationEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/LinearInterpolationEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/ArrayEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/QuadratureEvaluator.hpp>

#include <iostream>
#include <vector>

namespace sgpp {
namespace combigrid {

class CombigridMultiOperationImpl {
 public:
  CombigridMultiOperationImpl(
      std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>>> evaluatorPrototypes,
      std::shared_ptr<AbstractCombigridStorage> storage)
      : storage(storage),
        fullGridEval(new FullGridTensorEvaluator<FloatArrayVector>(storage, evaluatorPrototypes,
                                                                   pointHierarchies)),
        combiEval(new CombigridEvaluator<FloatArrayVector>(pointHierarchies.size(), fullGridEval)) {
  }

  std::shared_ptr<AbstractCombigridStorage> storage;
  std::shared_ptr<FullGridTensorEvaluator<FloatArrayVector>> fullGridEval;
  std::shared_ptr<CombigridEvaluator<FloatArrayVector>> combiEval;
};

CombigridMultiOperation::CombigridMultiOperation(
    std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
    std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>>> evaluatorPrototypes,
    MultiFunction func)
    : impl(new CombigridMultiOperationImpl(pointHierarchies, evaluatorPrototypes,
                                           std::shared_ptr<AbstractCombigridStorage>(
                                               new CombigridTreeStorage(pointHierarchies, func)))) {
}

CombigridMultiOperation::CombigridMultiOperation(
    std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
    std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>>> evaluatorPrototypes,
    std::shared_ptr<AbstractCombigridStorage> storage)
    : impl(new CombigridMultiOperationImpl(pointHierarchies, evaluatorPrototypes, storage)) {}

base::DataVector CombigridMultiOperation::evaluate(size_t q,
                                                   std::vector<base::DataVector> const &params) {
  // TODO(holzmudd): check dimensionalities and sizes...
  std::vector<FloatArrayVector> vecs(params[0].getSize());

  for (size_t i = 0; i < params.size(); ++i) {
    auto &param = params[i];
    for (size_t j = 0; j < param.getSize(); ++j) {
      vecs[j].at(i) = param[j];
    }
  }

  impl->fullGridEval->setParameters(vecs);
  impl->combiEval->clear();
  impl->combiEval->addRegularLevels(q);

  auto values = impl->combiEval->getValue().getValues();

  base::DataVector result(values.size());

  for (size_t i = 0; i < values.size(); ++i) {
    result[i] = values[i].getValue();
  }

  std::cout << "Storage Entries: " << impl->fullGridEval->getStorage()->getNumEntries() << "\n";

  return result;
}

base::DataVector CombigridMultiOperation::evaluateAdaptive(
    size_t maxNumPoints, std::vector<base::DataVector> const &params) {
  // TODO(holzmudd): check dimensionalities and sizes...
  std::vector<FloatArrayVector> vecs(params[0].getSize());

  for (size_t i = 0; i < params.size(); ++i) {
    auto &param = params[i];
    for (size_t j = 0; j < param.getSize(); ++j) {
      vecs[j].at(i) = param[j];
    }
  }

  impl->fullGridEval->setParameters(vecs);
  impl->combiEval->clear();
  impl->combiEval->addLevelsAdaptive(maxNumPoints);

  auto values = impl->combiEval->getValue().getValues();

  base::DataVector result(values.size());

  for (size_t i = 0; i < values.size(); ++i) {
    result[i] = values[i].getValue();
  }

  // std::cout << "Storage Entries: " << impl->fullGridEval->getStorage()->getNumEntries() << "\n";

  return result;
}

std::shared_ptr<AbstractMultiStorage<FloatArrayVector>> CombigridMultiOperation::getDifferences() {
  return impl->combiEval->differences();
}

std::shared_ptr<CombigridMultiOperation>
CombigridMultiOperation::createExpClenshawCurtisPolynomialInterpolation(size_t numDimensions,
                                                                        MultiFunction func) {
  return std::make_shared<CombigridMultiOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(
          numDimensions, std::make_shared<NestedPointHierarchy>(
                             std::make_shared<ClenshawCurtisDistribution>(),
                             std::make_shared<ExponentialLevelorderPointOrdering>())),
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>>>(
          numDimensions, std::make_shared<ArrayEvaluator<BarycentricInterpolationEvaluator>>(true)),
      func);
}

std::shared_ptr<CombigridMultiOperation>
CombigridMultiOperation::createLinearClenshawCurtisPolynomialInterpolation(size_t numDimensions,
                                                                           MultiFunction func) {
  return std::make_shared<CombigridMultiOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(
          numDimensions, std::make_shared<NonNestedPointHierarchy>(
                             std::make_shared<ClenshawCurtisDistribution>(),
                             std::make_shared<IdentityPointOrdering>(
                                 std::make_shared<LinearGrowthStrategy>(2), true))),
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>>>(
          numDimensions, std::make_shared<ArrayEvaluator<BarycentricInterpolationEvaluator>>(true)),
      func);
}

std::shared_ptr<CombigridMultiOperation>
CombigridMultiOperation::createExpLejaPolynomialInterpolation(size_t numDimensions,
                                                              MultiFunction func) {
  return std::make_shared<CombigridMultiOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(
          numDimensions, std::make_shared<NestedPointHierarchy>(
                             std::make_shared<LejaPointDistribution>(),
                             std::make_shared<IdentityPointOrdering>(
                                 std::make_shared<ExponentialGrowthStrategy>(), false))),
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>>>(
          numDimensions, std::make_shared<ArrayEvaluator<BarycentricInterpolationEvaluator>>(true)),
      func);
}

std::shared_ptr<CombigridMultiOperation>
CombigridMultiOperation::createLinearLejaPolynomialInterpolation(size_t numDimensions,
                                                                 MultiFunction func,
                                                                 size_t growthFactor) {
  return std::make_shared<CombigridMultiOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(
          numDimensions, std::make_shared<NestedPointHierarchy>(
                             std::make_shared<LejaPointDistribution>(),
                             std::make_shared<IdentityPointOrdering>(
                                 std::make_shared<LinearGrowthStrategy>(growthFactor), false))),
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>>>(
          numDimensions, std::make_shared<ArrayEvaluator<BarycentricInterpolationEvaluator>>(true)),
      func);
}

std::shared_ptr<CombigridMultiOperation> CombigridMultiOperation::createLinearLejaQuadrature(
    size_t numDimensions, MultiFunction func, size_t growthFactor) {
  return std::make_shared<CombigridMultiOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(
          numDimensions, std::make_shared<NestedPointHierarchy>(
                             std::make_shared<LejaPointDistribution>(),
                             std::make_shared<IdentityPointOrdering>(
                                 std::make_shared<LinearGrowthStrategy>(growthFactor), false))),
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>>>(
          numDimensions, std::make_shared<ArrayEvaluator<QuadratureEvaluator>>(false)),
      func);
}

std::shared_ptr<CombigridMultiOperation>
CombigridMultiOperation::createExpUniformLinearInterpolation(size_t numDimensions,
                                                             MultiFunction func) {
  return std::make_shared<CombigridMultiOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(
          numDimensions, std::make_shared<NestedPointHierarchy>(
                             std::make_shared<UniformPointDistribution>(),
                             std::make_shared<ExponentialLevelorderPointOrdering>())),
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>>>(
          numDimensions, std::make_shared<ArrayEvaluator<LinearInterpolationEvaluator>>(true)),
      func);
}

} /* namespace combigrid */
} /* namespace sgpp*/
