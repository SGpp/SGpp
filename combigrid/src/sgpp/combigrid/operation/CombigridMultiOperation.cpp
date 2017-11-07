// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/grid/distribution/ClenshawCurtisDistribution.hpp>
#include <sgpp/combigrid/grid/distribution/LejaPointDistribution.hpp>
#include <sgpp/combigrid/grid/distribution/UniformBoundaryPointDistribution.hpp>
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
#include <sgpp/combigrid/operation/multidim/WeightedRatioLevelManager.hpp>
#include <sgpp/combigrid/operation/multidim/fullgrid/AbstractFullGridEvaluationStrategy.hpp>
#include <sgpp/combigrid/operation/multidim/fullgrid/AbstractFullGridSummationStrategy.hpp>
#include <sgpp/combigrid/operation/multidim/fullgrid/FullGridCallbackEvaluator.hpp>
#include <sgpp/combigrid/operation/multidim/fullgrid/FullGridGridBasedEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/ArrayEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/BSplineInterpolationEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/BSplineRoutines.hpp>
#include <sgpp/combigrid/operation/onedim/LinearInterpolationEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/PolynomialInterpolationEvaluator.hpp>
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
      std::shared_ptr<LevelManager> levelManager, std::shared_ptr<AbstractCombigridStorage> storage,
      FullGridSummationStrategyType summationStrategyType = FullGridSummationStrategyType::LINEAR)
      : storage(storage), pointHierarchies(pointHierarchies), levelManager(levelManager) {
    fullGridEval = std::make_shared<FullGridCallbackEvaluator<FloatArrayVector>>(
        storage, evaluatorPrototypes, pointHierarchies, summationStrategyType);
    combiEval = std::make_shared<CombigridEvaluator<FloatArrayVector>>(pointHierarchies.size(),
                                                                       fullGridEval);

    levelManager->setLevelEvaluator(combiEval);
  }

  CombigridMultiOperationImpl(
      std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>>> evaluatorPrototypes,
      std::shared_ptr<LevelManager> levelManager, std::shared_ptr<AbstractCombigridStorage> storage,
      GridFunction gridFunc,
      FullGridSummationStrategyType summationStrategyType = FullGridSummationStrategyType::LINEAR)
      : storage(storage), pointHierarchies(pointHierarchies), levelManager(levelManager) {
    fullGridEval = std::make_shared<FullGridGridBasedEvaluator<FloatArrayVector>>(
        storage, evaluatorPrototypes, pointHierarchies, gridFunc, summationStrategyType);
    combiEval = std::make_shared<CombigridEvaluator<FloatArrayVector>>(pointHierarchies.size(),
                                                                       fullGridEval);
    levelManager->setLevelEvaluator(combiEval);
  }

  std::shared_ptr<AbstractCombigridStorage> storage;
  std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies;
  std::shared_ptr<AbstractFullGridEvaluationStrategy<FloatArrayVector>> fullGridEval;
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

CombigridMultiOperation::CombigridMultiOperation(
    std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
    std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>>> evaluatorPrototypes,
    std::shared_ptr<LevelManager> levelManager, GridFunction gridFunc, bool exploitNesting)
    : impl(new CombigridMultiOperationImpl(
          pointHierarchies, evaluatorPrototypes, levelManager,
          std::shared_ptr<AbstractCombigridStorage>(
              new CombigridTreeStorage(pointHierarchies, exploitNesting)),
          gridFunc)) {}

void CombigridMultiOperation::setParameters(const std::vector<base::DataVector> &params) {
  if (params.size() == 0) {
    throw std::runtime_error("CombigridMultiOperation::setParameters(): params.size() == 0");
  }
  size_t numParametersPerEvaluation = params[0].getSize();
  std::vector<FloatArrayVector> vecs(numParametersPerEvaluation);

  for (size_t i = 0; i < params.size(); ++i) {
    auto &param = params[i];

    if (param.getSize() != numParametersPerEvaluation) {
      throw std::runtime_error(
          "CombigridMultiOperation::setParameters(): not all parameter vectors have the same "
          "length");
    }
    for (size_t j = 0; j < numParametersPerEvaluation; ++j) {
      vecs[j].at(i) = param[j];
    }
  }

  impl->fullGridEval->setParameters(vecs);
  impl->combiEval->clear();
}

void CombigridMultiOperation::setParameters(const base::DataMatrix &params) {
  std::vector<FloatArrayVector> vecs(params.getNrows());

  for (size_t i = 0; i < params.getNrows(); ++i) {
    for (size_t j = 0; j < params.getNcols(); ++j) {
      vecs[i].at(j) = params(i, j);
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

std::shared_ptr<AbstractCombigridStorage> CombigridMultiOperation::getStorage() {
  return impl->storage;
}

std::vector<std::shared_ptr<AbstractPointHierarchy>>
CombigridMultiOperation::getPointHierarchies() {
  return impl->pointHierarchies;
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

  impl->levelManager->addRegularLevels(q);

  return getResult();
}

base::DataVector CombigridMultiOperation::evaluate(size_t q, base::DataMatrix const &params) {
  setParameters(params);

  impl->levelManager->addRegularLevels(q);

  return getResult();
}

std::shared_ptr<AbstractMultiStorage<FloatArrayVector>> CombigridMultiOperation::getDifferences() {
  return impl->combiEval->differences();
}

size_t CombigridMultiOperation::numStoredFunctionValues() { return impl->storage->getNumEntries(); }

size_t CombigridMultiOperation::numGridPoints() { return impl->levelManager->numGridPoints(); }

size_t CombigridMultiOperation::getUpperPointBound() const {
  return impl->levelManager->getUpperPointBound();
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
CombigridMultiOperation::createExpChebyshevPolynomialInterpolation(size_t numDimensions,
                                                                   MultiFunction func) {
  return std::make_shared<CombigridMultiOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(numDimensions,
                                                           CombiHierarchies::expChebyshev()),
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
CombigridMultiOperation::createExpL2LejaPolynomialInterpolation(size_t numDimensions,
                                                                MultiFunction func) {
  return std::make_shared<CombigridMultiOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(numDimensions,
                                                           CombiHierarchies::expL2Leja()),
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

std::shared_ptr<CombigridMultiOperation>
CombigridMultiOperation::createLinearL2LejaPolynomialInterpolation(size_t numDimensions,
                                                                   MultiFunction func,
                                                                   size_t growthFactor) {
  return std::make_shared<CombigridMultiOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(
          numDimensions, CombiHierarchies::linearL2Leja(growthFactor)),
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

std::shared_ptr<CombigridMultiOperation> CombigridMultiOperation::createLinearL2LejaQuadrature(
    size_t numDimensions, MultiFunction func, size_t growthFactor) {
  return std::make_shared<CombigridMultiOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(
          numDimensions, CombiHierarchies::linearL2Leja(growthFactor)),
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>>>(
          numDimensions, CombiEvaluators::multiQuadrature()),
      std::make_shared<StandardLevelManager>(), func);
}

std::shared_ptr<CombigridMultiOperation> CombigridMultiOperation::createExpClenshawCurtisQuadrature(
    size_t numDimensions, MultiFunction func) {
  return std::make_shared<CombigridMultiOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(numDimensions,
                                                           CombiHierarchies::expClenshawCurtis()),
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

std::shared_ptr<CombigridMultiOperation>
CombigridMultiOperation::createExpUniformBoundaryLinearInterpolation(size_t numDimensions,
                                                                     MultiFunction func) {
  return std::make_shared<CombigridMultiOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(numDimensions,
                                                           CombiHierarchies::expUniformBoundary()),
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>>>(
          numDimensions, CombiEvaluators::multiLinearInterpolation()),
      std::make_shared<StandardLevelManager>(), func);
}

std::shared_ptr<CombigridMultiOperation>
CombigridMultiOperation::createExpUniformBoundaryBsplineInterpolation(size_t numDimensions,
                                                                      MultiFunction func,
                                                                      size_t degree) {
  std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies(
      numDimensions, sgpp::combigrid::CombiHierarchies::expUniformBoundary());
  std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>>> evaluators(
      numDimensions, sgpp::combigrid::CombiEvaluators::multiBSplineInterpolation(degree));
  std::shared_ptr<sgpp::combigrid::LevelManager> levelManager(
      new sgpp::combigrid::WeightedRatioLevelManager());
  sgpp::combigrid::GridFunction gf = BSplineCoefficientGridFunction(func, pointHierarchies, degree);
  bool exploitNesting = false;
  return std::make_shared<sgpp::combigrid::CombigridMultiOperation>(
      pointHierarchies, evaluators, levelManager, gf, exploitNesting);
}

} /* namespace combigrid */
} /* namespace sgpp*/
