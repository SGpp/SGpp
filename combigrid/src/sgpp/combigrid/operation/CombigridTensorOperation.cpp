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
#include <sgpp/combigrid/operation/CombigridTensorOperation.hpp>
#include <sgpp/combigrid/operation/Configurations.hpp>
#include <sgpp/combigrid/operation/multidim/AveragingLevelManager.hpp>
#include <sgpp/combigrid/operation/multidim/CombigridEvaluator.hpp>
#include <sgpp/combigrid/operation/multidim/fullgrid/AbstractFullGridEvaluator.hpp>
#include <sgpp/combigrid/operation/multidim/fullgrid/AbstractFullGridSummationStrategy.hpp>
#include <sgpp/combigrid/operation/multidim/fullgrid/FullGridCallbackEvaluator.hpp>
#include <sgpp/combigrid/operation/multidim/fullgrid/FullGridGridBasedEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/ArrayEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/LinearInterpolationEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/PolynomialInterpolationEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/PolynomialQuadratureEvaluator.hpp>
#include <sgpp/combigrid/storage/tree/CombigridTreeStorage.hpp>
#include <sgpp/combigrid/operation/CombigridOperation.hpp>
#include <sgpp/combigrid/functions/OrthogonalPolynomialBasis1D.hpp>
#include <sgpp/combigrid/algebraic/NormStrategy.hpp>
#include <sgpp/combigrid/functions/OrthogonalBasisFunctionsCollection.hpp>
#include <sgpp/combigrid/utils/BSplineRoutines.hpp>

#include <sgpp/base/exception/generation_exception.hpp>

#include <iostream>
#include <vector>

typedef sgpp::combigrid::AveragingLevelManager StandardLevelManager;  // TODO(holzmudd)

namespace sgpp {
namespace combigrid {

class CombigridTensorOperationImpl {
 public:
  CombigridTensorOperationImpl(
      std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatTensorVector>>> evaluatorPrototypes,
      std::shared_ptr<LevelManager> levelManager, std::shared_ptr<AbstractCombigridStorage> storage,
      FullGridSummationStrategyType summationStrategyType,
      std::shared_ptr<NormStrategy<FloatTensorVector>> normStrategy)
      : storage(storage), pointHierarchies(pointHierarchies), levelManager(levelManager) {
    fullGridEval = std::make_shared<FullGridCallbackEvaluator<FloatTensorVector>>(
        storage, evaluatorPrototypes, pointHierarchies, summationStrategyType);
    combiEval = std::make_shared<CombigridEvaluator<FloatTensorVector>>(pointHierarchies.size(),
                                                                        fullGridEval, normStrategy);
    levelManager->setLevelEvaluator(combiEval);
  }

  CombigridTensorOperationImpl(
      std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatTensorVector>>> evaluatorPrototypes,
      std::shared_ptr<LevelManager> levelManager, std::shared_ptr<AbstractCombigridStorage> storage,
      GridFunction gridFunc, FullGridSummationStrategyType summationStrategyType,
      std::shared_ptr<NormStrategy<FloatTensorVector>> normStrategy)
      : storage(storage), pointHierarchies(pointHierarchies), levelManager(levelManager) {
    fullGridEval = std::make_shared<FullGridGridBasedEvaluator<FloatTensorVector>>(
        storage, evaluatorPrototypes, pointHierarchies, gridFunc, summationStrategyType);
    combiEval = std::make_shared<CombigridEvaluator<FloatTensorVector>>(pointHierarchies.size(),
                                                                        fullGridEval, normStrategy);
    levelManager->setLevelEvaluator(combiEval);
  }

  std::shared_ptr<AbstractCombigridStorage> storage;
  std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies;
  std::shared_ptr<AbstractFullGridEvaluator<FloatTensorVector>> fullGridEval;
  std::shared_ptr<CombigridEvaluator<FloatTensorVector>> combiEval;
  std::shared_ptr<LevelManager> levelManager;
};

CombigridTensorOperation::CombigridTensorOperation(
    std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
    std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatTensorVector>>> evaluatorPrototypes,
    std::shared_ptr<LevelManager> levelManager, MultiFunction func, bool exploitNesting,
    FullGridSummationStrategyType summationStrategyType,
    std::shared_ptr<NormStrategy<FloatTensorVector>> normStrategy) {
  impl = std::make_shared<CombigridTensorOperationImpl>(
      pointHierarchies, evaluatorPrototypes, levelManager,
      std::make_shared<CombigridTreeStorage>(pointHierarchies, exploitNesting, func),
      summationStrategyType, normStrategy);
}

CombigridTensorOperation::CombigridTensorOperation(
    std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
    std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatTensorVector>>> evaluatorPrototypes,
    std::shared_ptr<LevelManager> levelManager, std::shared_ptr<AbstractCombigridStorage> storage,
    FullGridSummationStrategyType summationStrategyType,
    std::shared_ptr<NormStrategy<FloatTensorVector>> normStrategy) {
  impl = std::make_shared<CombigridTensorOperationImpl>(pointHierarchies, evaluatorPrototypes,
                                                        levelManager, storage,
                                                        summationStrategyType, normStrategy);
}

CombigridTensorOperation::CombigridTensorOperation(
    std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
    std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatTensorVector>>> evaluatorPrototypes,
    std::shared_ptr<LevelManager> levelManager, GridFunction gridFunc, bool exploitNesting,
    FullGridSummationStrategyType summationStrategyType,
    std::shared_ptr<NormStrategy<FloatTensorVector>> normStrategy) {
  impl = std::make_shared<CombigridTensorOperationImpl>(
      pointHierarchies, evaluatorPrototypes, levelManager,
      std::make_shared<CombigridTreeStorage>(pointHierarchies, exploitNesting), gridFunc,
      summationStrategyType, normStrategy);
}

void CombigridTensorOperation::setParameters(const std::vector<FloatTensorVector> &params) {
  impl->fullGridEval->setParameters(params);
  impl->combiEval->clear();
}

FloatTensorVector CombigridTensorOperation::getResult() { return impl->combiEval->getValue(); }

std::shared_ptr<AbstractCombigridStorage> CombigridTensorOperation::getStorage() {
  return impl->storage;
}

std::shared_ptr<LevelManager> CombigridTensorOperation::getLevelManager() {
  return impl->levelManager;
}

void CombigridTensorOperation::setLevelManager(std::shared_ptr<LevelManager> levelManager) {
  levelManager->setLevelEvaluator(impl->combiEval);
  impl->levelManager = levelManager;
}

FloatTensorVector CombigridTensorOperation::evaluate(size_t q,
                                                     std::vector<FloatTensorVector> const &params) {
  setParameters(params);

  impl->levelManager->addRegularLevels(q);

  return getResult();
}

std::shared_ptr<AbstractMultiStorage<FloatTensorVector>>
CombigridTensorOperation::getDifferences() {
  return impl->combiEval->differences();
}

size_t CombigridTensorOperation::numStoredFunctionValues() {
  return impl->storage->getNumEntries();
}

size_t CombigridTensorOperation::numGridPoints() { return impl->levelManager->numGridPoints(); }

size_t CombigridTensorOperation::numDims() { return impl->levelManager->numDims(); }

size_t CombigridTensorOperation::getUpperPointBound() const {
  return impl->levelManager->getUpperPointBound();
}

std::vector<std::shared_ptr<AbstractPointHierarchy>>
CombigridTensorOperation::getPointHierarchies() {
  return impl->pointHierarchies;
}

std::shared_ptr<CombigridTensorOperation>
CombigridTensorOperation::createExpClenshawCurtisPolynomialInterpolation(
    std::shared_ptr<AbstractInfiniteFunctionBasis1D> functionBasis, size_t numDimensions,
    MultiFunction func) {
  return std::make_shared<CombigridTensorOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(numDimensions,
                                                           CombiHierarchies::expClenshawCurtis()),
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatTensorVector>>>(
          numDimensions, CombiEvaluators::tensorInterpolation(functionBasis)),
      std::make_shared<StandardLevelManager>(), func);
}

std::shared_ptr<CombigridTensorOperation>
CombigridTensorOperation::createExpChebyshevPolynomialInterpolation(
    std::shared_ptr<AbstractInfiniteFunctionBasis1D> functionBasis, size_t numDimensions,
    MultiFunction func) {
  return std::make_shared<CombigridTensorOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(numDimensions,
                                                           CombiHierarchies::expChebyshev()),
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatTensorVector>>>(
          numDimensions, CombiEvaluators::tensorInterpolation(functionBasis)),
      std::make_shared<StandardLevelManager>(), func);
}

std::shared_ptr<CombigridTensorOperation>
CombigridTensorOperation::createLinearClenshawCurtisPolynomialInterpolation(
    std::shared_ptr<AbstractInfiniteFunctionBasis1D> functionBasis, size_t numDimensions,
    MultiFunction func) {
  return std::make_shared<CombigridTensorOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(
          numDimensions, std::make_shared<NonNestedPointHierarchy>(
                             std::make_shared<ClenshawCurtisDistribution>(),
                             std::make_shared<IdentityPointOrdering>(
                                 std::make_shared<LinearGrowthStrategy>(2), true))),
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatTensorVector>>>(
          numDimensions, CombiEvaluators::tensorInterpolation(functionBasis)),
      std::make_shared<StandardLevelManager>(), func);
}

std::shared_ptr<CombigridTensorOperation>
CombigridTensorOperation::createExpLejaPolynomialInterpolation(
    std::shared_ptr<AbstractInfiniteFunctionBasis1D> functionBasis, size_t numDimensions,
    MultiFunction func) {
  return std::make_shared<CombigridTensorOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(numDimensions,
                                                           CombiHierarchies::expLeja()),
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatTensorVector>>>(
          numDimensions, CombiEvaluators::tensorInterpolation(functionBasis)),
      std::make_shared<StandardLevelManager>(), func);
}

std::shared_ptr<CombigridTensorOperation>
CombigridTensorOperation::createExpL2LejaPolynomialInterpolation(
    std::shared_ptr<AbstractInfiniteFunctionBasis1D> functionBasis, size_t numDimensions,
    MultiFunction func) {
  return std::make_shared<CombigridTensorOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(numDimensions,
                                                           CombiHierarchies::expL2Leja()),
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatTensorVector>>>(
          numDimensions, CombiEvaluators::tensorInterpolation(functionBasis)),
      std::make_shared<StandardLevelManager>(), func);
}

std::shared_ptr<CombigridTensorOperation>
CombigridTensorOperation::createLinearLejaPolynomialInterpolation(
    std::shared_ptr<AbstractInfiniteFunctionBasis1D> functionBasis, size_t numDimensions,
    MultiFunction func, size_t growthFactor) {
  return std::make_shared<CombigridTensorOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(
          numDimensions, CombiHierarchies::linearLeja(growthFactor)),
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatTensorVector>>>(
          numDimensions, CombiEvaluators::tensorInterpolation(functionBasis)),
      std::make_shared<StandardLevelManager>(), func);
}

std::shared_ptr<CombigridTensorOperation>
CombigridTensorOperation::createLinearL2LejaPolynomialInterpolation(
    std::shared_ptr<AbstractInfiniteFunctionBasis1D> functionBasis, size_t numDimensions,
    MultiFunction func, size_t growthFactor) {
  return std::make_shared<CombigridTensorOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(
          numDimensions, CombiHierarchies::linearL2Leja(growthFactor)),
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatTensorVector>>>(
          numDimensions, CombiEvaluators::tensorInterpolation(functionBasis)),
      std::make_shared<StandardLevelManager>(), func);
}

std::shared_ptr<CombigridTensorOperation>
CombigridTensorOperation::createOperationTensorPolynomialInterpolation(
    std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
    std::shared_ptr<AbstractCombigridStorage> storage,
    std::shared_ptr<AbstractInfiniteFunctionBasis1D> functionBasis,
    FullGridSummationStrategyType summationStrategyType) {
  size_t numDims = pointHierarchies.size();
  std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatTensorVector>>> evaluatorPrototypes(
      numDims, sgpp::combigrid::CombiEvaluators::tensorInterpolation(functionBasis));

  auto tensorOperation = std::make_shared<CombigridTensorOperation>(
      pointHierarchies, evaluatorPrototypes, std::make_shared<StandardLevelManager>(), storage,
      summationStrategyType);

  return tensorOperation;
}

std::shared_ptr<CombigridTensorOperation>
CombigridTensorOperation::createOperationTensorPolynomialInterpolation(
    std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
    std::shared_ptr<AbstractCombigridStorage> storage,
    OrthogonalBasisFunctionsCollection &functionBases,
    FullGridSummationStrategyType summationStrategyType) {
  // check if dimensionalities match
  if (functionBases.size() != pointHierarchies.size()) {
    throw sgpp::base::generation_exception(
        "CombigridTensorOperation::createOperationTensorPolynomialInterpolation - the "
        "dimensionalities of the basis functions and the point hierarchies do not match");
  }

  // create combi evaluators
  sgpp::combigrid::CombiEvaluators::TensorCollection evaluators;
  for (size_t i = 0; i < functionBases.size(); i++) {
    evaluators.push_back(sgpp::combigrid::CombiEvaluators::tensorInterpolation(functionBases[i]));
  }

  auto tensorOperation = std::make_shared<CombigridTensorOperation>(
      pointHierarchies, evaluators, std::make_shared<StandardLevelManager>(), storage,
      summationStrategyType);

  return tensorOperation;
}

std::shared_ptr<CombigridTensorOperation>
CombigridTensorOperation::createOperationTensorBSplineInterpolation(
    std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
    std::shared_ptr<AbstractCombigridStorage> storage, std::shared_ptr<LevelManager> levelManager,
    FullGridSummationStrategyType summationStrategyType) {
  // create combi evaluators
  size_t numDims = pointHierarchies.size();
  sgpp::combigrid::CombiEvaluators::TensorCollection evaluators;
  for (size_t i = 0; i < numDims; i++) {
    evaluators.push_back(sgpp::combigrid::CombiEvaluators::tensorBSplineInterpolation());
  }

  auto tensorOperation = std::make_shared<CombigridTensorOperation>(
      pointHierarchies, evaluators, std::make_shared<StandardLevelManager>(), storage,
      summationStrategyType);
  auto levelStructure = levelManager->getLevelStructure();
  tensorOperation->getLevelManager()->addLevelsFromStructure(levelStructure);

  return tensorOperation;
}

std::shared_ptr<CombigridTensorOperation>
CombigridTensorOperation::createExpUniformBoundaryBSplineInterpolation(size_t numDimensions,
                                                                       MultiFunction func,
                                                                       size_t degree) {
  sgpp::combigrid::CombiHierarchies::Collection pointHierarchies(
      numDimensions, sgpp::combigrid::CombiHierarchies::expUniformBoundary());
  sgpp::combigrid::CombiEvaluators::TensorCollection evaluators(
      numDimensions, sgpp::combigrid::CombiEvaluators::tensorBSplineInterpolation(degree));

  std::shared_ptr<sgpp::combigrid::LevelManager> levelManager =
      std::make_shared<sgpp::combigrid::AveragingLevelManager>();
  sgpp::combigrid::GridFunction gf =
      BSplineTensorCoefficientGridFunction(func, pointHierarchies, degree);
  bool exploitNesting = false;
  return std::make_shared<sgpp::combigrid::CombigridTensorOperation>(
      pointHierarchies, evaluators, levelManager, gf, exploitNesting);
}

} /* namespace combigrid */
} /* namespace sgpp*/
