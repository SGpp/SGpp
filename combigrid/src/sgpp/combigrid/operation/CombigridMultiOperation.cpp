// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/grid/distribution/ClenshawCurtisDistribution.hpp>
#include <sgpp/combigrid/grid/distribution/LejaPointDistribution.hpp>
#include <sgpp/combigrid/grid/distribution/UniformPointDistribution.hpp>
#include <sgpp/combigrid/grid/growth/ExponentialGrowthStrategy.hpp>
#include <sgpp/combigrid/grid/growth/LinearGrowthStrategy.hpp>
#include <sgpp/combigrid/grid/hierarchy/NestedPointHierarchy.hpp>
#include <sgpp/combigrid/grid/hierarchy/NonNestedPointHierarchy.hpp>
#include <sgpp/combigrid/grid/ordering/ExponentialLevelorderPointOrdering.hpp>
#include <sgpp/combigrid/grid/ordering/IdentityPointOrdering.hpp>
#include <sgpp/combigrid/operation/CombigridMultiOperation.hpp>
#include <sgpp/combigrid/operation/Configurations.hpp>
#include <sgpp/combigrid/operation/multidim/AveragingLevelManager.hpp>
#include <sgpp/combigrid/operation/multidim/CombigridEvaluator.hpp>
#include <sgpp/combigrid/operation/multidim/FullGridTensorEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/ArrayEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/BarycentricInterpolationEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/LinearInterpolationEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/QuadratureEvaluator.hpp>
#include <sgpp/combigrid/storage/tree/CombigridTreeStorage.hpp>

#include <iostream>
#include <vector>

typedef sgpp::combigrid::AveragingLevelManager StandardLevelManager;  // TODO(holzmudd)

namespace sgpp {
namespace combigrid {

class CombigridMultiOperationImpl {
 public:
  CombigridMultiOperationImpl(
      std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>>> evaluatorPrototypes,
      std::shared_ptr<LevelManager> levelManager, std::shared_ptr<AbstractCombigridStorage> storage)
      : storage(storage),
        fullGridEval(new FullGridTensorEvaluator<FloatArrayVector>(storage, evaluatorPrototypes,
                                                                   pointHierarchies)),
        combiEval(new CombigridEvaluator<FloatArrayVector>(pointHierarchies.size(), fullGridEval)),
        levelManager(levelManager) {
    levelManager->setLevelEvaluator(combiEval);
  }

  std::shared_ptr<AbstractCombigridStorage> storage;
  std::shared_ptr<FullGridTensorEvaluator<FloatArrayVector>> fullGridEval;
  std::shared_ptr<CombigridEvaluator<FloatArrayVector>> combiEval;
  std::shared_ptr<LevelManager> levelManager;
};

CombigridMultiOperation::CombigridMultiOperation(
    std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
    std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>>> evaluatorPrototypes,
    std::shared_ptr<LevelManager> levelManager, MultiFunction func)
    : impl(new CombigridMultiOperationImpl(pointHierarchies, evaluatorPrototypes, levelManager,
                                           std::shared_ptr<AbstractCombigridStorage>(
                                               new CombigridTreeStorage(pointHierarchies, func)))) {
}

CombigridMultiOperation::CombigridMultiOperation(
    std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
    std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>>> evaluatorPrototypes,
    std::shared_ptr<LevelManager> levelManager, std::shared_ptr<AbstractCombigridStorage> storage)
    : impl(new CombigridMultiOperationImpl(pointHierarchies, evaluatorPrototypes, levelManager,
                                           storage)) {}

void CombigridMultiOperation::setParameters(const std::vector<base::DataVector> &params) {
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
}

void CombigridMultiOperation::setParameters(const base::DataMatrix &params) {
  std::vector<FloatArrayVector> vecs(params.getNcols());

  for (size_t i = 0; i < params.getNrows(); ++i) {
    for (size_t j = 0; j < params.getNcols(); ++j) {
      vecs[j].at(i) = params(i, j);
    }
  }

  impl->fullGridEval->setParameters(vecs);
  impl->combiEval->clear();
}

base::DataVector CombigridMultiOperation::getResult() {
  auto values = impl->combiEval->getValue().getValues();

  base::DataVector result(values.size());

  for (size_t i = 0; i < values.size(); ++i) {
    result[i] = values[i].getValue();
  }

  return result;
}

std::shared_ptr<LevelManager> CombigridMultiOperation::getLevelManager() {
  return impl->levelManager;
}

void CombigridMultiOperation::setLevelManager(std::shared_ptr<LevelManager> levelManager) {
  levelManager->setLevelEvaluator(impl->combiEval);
  impl->levelManager = levelManager;
}

base::DataVector CombigridMultiOperation::evaluate(size_t q,
                                                   std::vector<base::DataVector> const &params) {
  setParameters(params);

  impl->combiEval->addRegularLevels(q);

  return getResult();
}

std::shared_ptr<AbstractMultiStorage<FloatArrayVector>> CombigridMultiOperation::getDifferences() {
  return impl->combiEval->differences();
}

std::shared_ptr<CombigridMultiOperation>
CombigridMultiOperation::createExpClenshawCurtisPolynomialInterpolation(size_t numDimensions,
                                                                        MultiFunction func) {
  return std::make_shared<CombigridMultiOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(numDimensions,
                                                           CombiHierarchies::expClenshawCurtis()),
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>>>(
          numDimensions, CombiEvaluators::multiPolynomialInterpolation()),
      std::make_shared<StandardLevelManager>(), func);
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
          numDimensions, CombiEvaluators::multiPolynomialInterpolation()),
      std::make_shared<StandardLevelManager>(), func);
}

std::shared_ptr<CombigridMultiOperation>
CombigridMultiOperation::createExpLejaPolynomialInterpolation(size_t numDimensions,
                                                              MultiFunction func) {
  return std::make_shared<CombigridMultiOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(numDimensions,
                                                           CombiHierarchies::expLeja()),
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>>>(
          numDimensions, CombiEvaluators::multiPolynomialInterpolation()),
      std::make_shared<StandardLevelManager>(), func);
}

std::shared_ptr<CombigridMultiOperation>
CombigridMultiOperation::createLinearLejaPolynomialInterpolation(size_t numDimensions,
                                                                 MultiFunction func,
                                                                 size_t growthFactor) {
  return std::make_shared<CombigridMultiOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(
          numDimensions, CombiHierarchies::linearLeja(growthFactor)),
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>>>(
          numDimensions, CombiEvaluators::multiPolynomialInterpolation()),
      std::make_shared<StandardLevelManager>(), func);
}

std::shared_ptr<CombigridMultiOperation> CombigridMultiOperation::createLinearLejaQuadrature(
    size_t numDimensions, MultiFunction func, size_t growthFactor) {
  return std::make_shared<CombigridMultiOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(
          numDimensions, CombiHierarchies::linearLeja(growthFactor)),
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>>>(
          numDimensions, CombiEvaluators::multiQuadrature()),
      std::make_shared<StandardLevelManager>(), func);
}

std::shared_ptr<CombigridMultiOperation>
CombigridMultiOperation::createExpUniformLinearInterpolation(size_t numDimensions,
                                                             MultiFunction func) {
  return std::make_shared<CombigridMultiOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(numDimensions,
                                                           CombiHierarchies::expUniform()),
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>>>(
          numDimensions, CombiEvaluators::multiLinearInterpolation()),
      std::make_shared<StandardLevelManager>(), func);
}

} /* namespace combigrid */
} /* namespace sgpp*/
