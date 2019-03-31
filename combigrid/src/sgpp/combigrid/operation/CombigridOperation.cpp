// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/grid/distribution/ClenshawCurtisDistribution.hpp>
#include <sgpp/combigrid/grid/distribution/UniformPointDistribution.hpp>
#include <sgpp/combigrid/grid/growth/ExponentialGrowthStrategy.hpp>
#include <sgpp/combigrid/grid/hierarchy/NestedPointHierarchy.hpp>
#include <sgpp/combigrid/grid/ordering/ExponentialLevelorderPointOrdering.hpp>
#include <sgpp/combigrid/operation/CombigridOperation.hpp>
#include <sgpp/combigrid/operation/multidim/AveragingLevelManager.hpp>
#include <sgpp/combigrid/operation/multidim/CombigridEvaluator.hpp>
#include <sgpp/combigrid/operation/multidim/WeightedRatioLevelManager.hpp>
#include <sgpp/combigrid/operation/multidim/fullgrid/AbstractFullGridEvaluationStrategy.hpp>
#include <sgpp/combigrid/operation/multidim/fullgrid/AbstractFullGridEvaluator.hpp>
#include <sgpp/combigrid/operation/multidim/fullgrid/AbstractFullGridSummationStrategy.hpp>
#include <sgpp/combigrid/operation/multidim/fullgrid/FullGridCallbackEvaluator.hpp>
#include <sgpp/combigrid/operation/multidim/fullgrid/FullGridGridBasedEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/LinearInterpolationEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/PolynomialInterpolationEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/PolynomialQuadratureEvaluator.hpp>
#include <sgpp/combigrid/storage/tree/CombigridTreeStorage.hpp>
#include <sgpp/combigrid/utils/BSplineRoutines.hpp>
#include <sgpp/combigrid/utils/CombigridBSplineBasis.hpp>
#include <sgpp/optimization/sle/solver/Auto.hpp>
#include <sgpp/optimization/sle/system/FullSLE.hpp>
#include <sgpp/optimization/tools/Printer.hpp>

#include <vector>

typedef sgpp::combigrid::AveragingLevelManager StandardLevelManager;

namespace sgpp {
namespace combigrid {

class CombigridOperationImpl {
 public:
  CombigridOperationImpl(
      std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>> evaluatorPrototypes,
      std::shared_ptr<LevelManager> levelManager, std::shared_ptr<AbstractCombigridStorage> storage,
      FullGridSummationStrategyType summationStrategyType,
      std::shared_ptr<NormStrategy<FloatScalarVector>> normStrategy)
      : storage(storage), pointHierarchies(pointHierarchies), levelManager(levelManager) {
    fullGridEval = std::make_shared<FullGridCallbackEvaluator<FloatScalarVector>>(
        storage, evaluatorPrototypes, pointHierarchies, summationStrategyType);
    combiEval = std::make_shared<CombigridEvaluator<FloatScalarVector>>(pointHierarchies.size(),
                                                                        fullGridEval, normStrategy);
    levelManager->setLevelEvaluator(combiEval);
  }

  CombigridOperationImpl(
      std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>> evaluatorPrototypes,
      std::shared_ptr<LevelManager> levelManager, std::shared_ptr<AbstractCombigridStorage> storage,
      GridFunction gridFunc, FullGridSummationStrategyType summationStrategyType,
      std::shared_ptr<NormStrategy<FloatScalarVector>> normStrategy)
      : storage(storage), pointHierarchies(pointHierarchies), levelManager(levelManager) {
    fullGridEval = std::make_shared<FullGridGridBasedEvaluator<FloatScalarVector>>(
        storage, evaluatorPrototypes, pointHierarchies, gridFunc, summationStrategyType);
    combiEval = std::make_shared<CombigridEvaluator<FloatScalarVector>>(pointHierarchies.size(),
                                                                        fullGridEval, normStrategy);
    levelManager->setLevelEvaluator(combiEval);
  }

  std::shared_ptr<AbstractCombigridStorage> storage;
  std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies;
  std::shared_ptr<AbstractFullGridEvaluator<FloatScalarVector>> fullGridEval;
  std::shared_ptr<CombigridEvaluator<FloatScalarVector>> combiEval;
  std::shared_ptr<LevelManager> levelManager;
};

CombigridOperation::CombigridOperation(
    std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
    std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>> evaluatorPrototypes,
    std::shared_ptr<LevelManager> levelManager, MultiFunction func, bool exploitNesting,
    FullGridSummationStrategyType summationStrategyType,
    std::shared_ptr<NormStrategy<FloatScalarVector>> normStrategy) {
  impl = std::make_shared<CombigridOperationImpl>(
      pointHierarchies, evaluatorPrototypes, levelManager,
      std::make_shared<CombigridTreeStorage>(pointHierarchies, exploitNesting, func),
      summationStrategyType, normStrategy);
}

CombigridOperation::CombigridOperation(
    std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
    std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>> evaluatorPrototypes,
    std::shared_ptr<LevelManager> levelManager, std::shared_ptr<AbstractCombigridStorage> storage,
    FullGridSummationStrategyType summationStrategyType,
    std::shared_ptr<NormStrategy<FloatScalarVector>> normStrategy) {
  impl =
      std::make_shared<CombigridOperationImpl>(pointHierarchies, evaluatorPrototypes, levelManager,
                                               storage, summationStrategyType, normStrategy);
}

CombigridOperation::CombigridOperation(
    std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
    std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>> evaluatorPrototypes,
    std::shared_ptr<LevelManager> levelManager, GridFunction gridFunc, bool exploitNesting,
    FullGridSummationStrategyType summationStrategyType,
    std::shared_ptr<NormStrategy<FloatScalarVector>> normStrategy) {
  impl = std::make_shared<CombigridOperationImpl>(
      pointHierarchies, evaluatorPrototypes, levelManager,
      std::make_shared<CombigridTreeStorage>(pointHierarchies, exploitNesting), gridFunc,
      summationStrategyType, normStrategy);
}

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
  impl->levelManager->addRegularLevels(q);

  return getResult();
}

double CombigridOperation::evaluateParallel(size_t q, base::DataVector const& param,
        size_t numThreads) {
  setParameters(param);
  impl->levelManager->addRegularLevelsParallel(q, numThreads);

  return getResult();
}

std::shared_ptr<AbstractCombigridStorage> CombigridOperation::getStorage() { return impl->storage; }

void CombigridOperation::setStorage(std::shared_ptr<AbstractCombigridStorage> storage) {
  impl->storage = storage;
}

std::vector<std::shared_ptr<AbstractPointHierarchy>> CombigridOperation::getPointHierarchies() {
  return impl->pointHierarchies;
}

std::shared_ptr<AbstractFullGridEvaluator<FloatScalarVector>>
CombigridOperation::getFullGridEval() {
  return impl->fullGridEval;
}

std::shared_ptr<LevelManager> CombigridOperation::getLevelManager() { return impl->levelManager; }

void CombigridOperation::setLevelManager(std::shared_ptr<LevelManager> levelManager) {
  levelManager->setLevelEvaluator(impl->combiEval);
  impl->levelManager = levelManager;
}

size_t CombigridOperation::numStoredFunctionValues() { return impl->storage->getNumEntries(); }

size_t CombigridOperation::numGridPoints() { return impl->levelManager->numGridPoints(); }

size_t CombigridOperation::numDims() { return impl->levelManager->numDims(); }

size_t CombigridOperation::getUpperPointBound() const {
  return impl->levelManager->getUpperPointBound();
}

std::shared_ptr<CombigridOperation>
CombigridOperation::createExpClenshawCurtisPolynomialInterpolation(size_t numDimensions,
                                                                   MultiFunction func) {
  return std::make_shared<CombigridOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(numDimensions,
                                                           CombiHierarchies::expClenshawCurtis()),
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>>(
          numDimensions, CombiEvaluators::polynomialInterpolation()),
      std::make_shared<StandardLevelManager>(), func);
}

std::shared_ptr<CombigridOperation> CombigridOperation::createExpChebyshevPolynomialInterpolation(
    size_t numDimensions, MultiFunction func) {
  return std::make_shared<CombigridOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(numDimensions,
                                                           CombiHierarchies::expChebyshev()),
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>>(
          numDimensions, CombiEvaluators::polynomialInterpolation()),
      std::make_shared<StandardLevelManager>(), func);
}

std::shared_ptr<CombigridOperation> CombigridOperation::createExpLejaPolynomialInterpolation(
    size_t numDimensions, MultiFunction func) {
  return std::make_shared<CombigridOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(numDimensions,
                                                           CombiHierarchies::expLeja()),
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>>(
          numDimensions, CombiEvaluators::polynomialInterpolation()),
      std::make_shared<StandardLevelManager>(), func);
}

std::shared_ptr<CombigridOperation> CombigridOperation::createExpL2LejaPolynomialInterpolation(
    size_t numDimensions, MultiFunction func) {
  return std::make_shared<CombigridOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(numDimensions,
                                                           CombiHierarchies::expL2Leja()),
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>>(
          numDimensions, CombiEvaluators::polynomialInterpolation()),
      std::make_shared<StandardLevelManager>(), func);
}

std::shared_ptr<CombigridOperation> CombigridOperation::createExpUniformPolynomialInterpolation(
    size_t numDimensions, MultiFunction func) {
  return std::make_shared<CombigridOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(numDimensions,
                                                           CombiHierarchies::expUniform()),
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>>(
          numDimensions, CombiEvaluators::polynomialInterpolation()),
      std::make_shared<StandardLevelManager>(), func);
}

std::shared_ptr<CombigridOperation>
CombigridOperation::createExpUniformBoundaryPolynomialInterpolation(size_t numDimensions,
                                                                    MultiFunction func) {
  return std::make_shared<CombigridOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(numDimensions,
                                                           CombiHierarchies::expUniformBoundary()),
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>>(
          numDimensions, CombiEvaluators::polynomialInterpolation()),
      std::make_shared<StandardLevelManager>(), func);
}

std::shared_ptr<CombigridOperation> CombigridOperation::createExpUniformLinearInterpolation(
    size_t numDimensions, MultiFunction func) {
  return std::make_shared<CombigridOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(numDimensions,
                                                           CombiHierarchies::expUniform()),
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>>(
          numDimensions, CombiEvaluators::linearInterpolation()),
      std::make_shared<StandardLevelManager>(), func);
}

std::shared_ptr<CombigridOperation> CombigridOperation::createExpUniformBoundaryLinearInterpolation(
    size_t numDimensions, MultiFunction func) {
  return std::make_shared<CombigridOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(numDimensions,
                                                           CombiHierarchies::expUniformBoundary()),
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>>(
          numDimensions, CombiEvaluators::linearInterpolation()),
      std::make_shared<StandardLevelManager>(), func);
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
      std::make_shared<StandardLevelManager>(), func);
}

std::shared_ptr<CombigridOperation> CombigridOperation::createLinearLejaPolynomialInterpolation(
    size_t numDimensions, MultiFunction func, size_t growthFactor) {
  return std::make_shared<CombigridOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(
          numDimensions, CombiHierarchies::linearLeja(growthFactor)),
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>>(
          numDimensions, CombiEvaluators::polynomialInterpolation()),
      std::make_shared<StandardLevelManager>(), func);
}

std::shared_ptr<CombigridOperation> CombigridOperation::createLinearL2LejaPolynomialInterpolation(
    size_t numDimensions, MultiFunction func, size_t growthFactor) {
  return std::make_shared<CombigridOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(
          numDimensions, CombiHierarchies::linearL2Leja(growthFactor)),
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>>(
          numDimensions, CombiEvaluators::polynomialInterpolation()),
      std::make_shared<StandardLevelManager>(), func);
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
      std::make_shared<StandardLevelManager>(), func);
}

std::shared_ptr<CombigridOperation>
CombigridOperation::createLinearUniformBoundaryPolynomialInterpolation(size_t numDimensions,
                                                                       MultiFunction func) {
  return std::make_shared<CombigridOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(
          numDimensions, std::make_shared<NonNestedPointHierarchy>(
                             std::make_shared<UniformBoundaryPointDistribution>(),
                             std::make_shared<IdentityPointOrdering>(
                                 std::make_shared<LinearGrowthStrategy>(2), true))),
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>>(
          numDimensions, CombiEvaluators::polynomialInterpolation()),
      std::make_shared<StandardLevelManager>(), func);
}

std::shared_ptr<CombigridOperation> CombigridOperation::createLinearLejaQuadrature(
    size_t numDimensions, MultiFunction func, size_t growthFactor) {
  return std::make_shared<CombigridOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(
          numDimensions, CombiHierarchies::linearLeja(growthFactor)),
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>>(
          numDimensions, CombiEvaluators::quadrature()),
      std::make_shared<StandardLevelManager>(), func);
}

std::shared_ptr<CombigridOperation> CombigridOperation::createLinearL2LejaQuadrature(
    size_t numDimensions, MultiFunction func, size_t growthFactor) {
  return std::make_shared<CombigridOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(
          numDimensions, CombiHierarchies::linearL2Leja(growthFactor)),
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>>(
          numDimensions, CombiEvaluators::quadrature()),
      std::make_shared<StandardLevelManager>(), func);
}

std::shared_ptr<CombigridOperation> CombigridOperation::createExpClenshawCurtisQuadrature(
    size_t numDimensions, MultiFunction func) {
  return std::make_shared<CombigridOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(numDimensions,
                                                           CombiHierarchies::expClenshawCurtis()),
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>>(
          numDimensions, CombiEvaluators::quadrature()),
      std::make_shared<StandardLevelManager>(), func);
}

std::shared_ptr<CombigridOperation> CombigridOperation::auxiliaryBsplineFunction(
    size_t numDimensions, MultiFunction func,
    sgpp::combigrid::CombiHierarchies::Collection pointHierarchies,
    sgpp::combigrid::CombiEvaluators::Collection evaluators, size_t degree) {
  // So far only WeightedRatioLevelManager has been used
  std::shared_ptr<sgpp::combigrid::LevelManager> levelManager =
      std::make_shared<sgpp::combigrid::AveragingLevelManager>();
  sgpp::combigrid::GridFunction gf = BSplineCoefficientGridFunction(func, pointHierarchies, degree);
  bool exploitNesting = false;
  return std::make_shared<sgpp::combigrid::CombigridOperation>(pointHierarchies, evaluators,
                                                               levelManager, gf, exploitNesting);
}

std::shared_ptr<CombigridOperation>
CombigridOperation::createExpUniformBoundaryBsplineInterpolation(size_t numDimensions,
                                                                 MultiFunction func,
                                                                 size_t degree = 3) {
  sgpp::combigrid::CombiHierarchies::Collection pointHierarchies(
      numDimensions, sgpp::combigrid::CombiHierarchies::expUniformBoundary());
  sgpp::combigrid::CombiEvaluators::Collection evaluators(
      numDimensions, sgpp::combigrid::CombiEvaluators::BSplineInterpolation(degree));
  return auxiliaryBsplineFunction(numDimensions, func, pointHierarchies, evaluators, degree);
}

std::shared_ptr<CombigridOperation> CombigridOperation::createExpClenshawCurtisBsplineInterpolation(
    size_t numDimensions, MultiFunction func, size_t degree = 3) {
  sgpp::combigrid::CombiHierarchies::Collection pointHierarchies(
      numDimensions, sgpp::combigrid::CombiHierarchies::expClenshawCurtis());
  sgpp::combigrid::CombiEvaluators::Collection evaluators(
      numDimensions, sgpp::combigrid::CombiEvaluators::BSplineInterpolation(degree));
  return auxiliaryBsplineFunction(numDimensions, func, pointHierarchies, evaluators, degree);
}

std::shared_ptr<CombigridOperation> CombigridOperation::createExpChebyshevBsplineInterpolation(
    size_t numDimensions, MultiFunction func, size_t degree = 3) {
  sgpp::combigrid::CombiHierarchies::Collection pointHierarchies(
      numDimensions, sgpp::combigrid::CombiHierarchies::expChebyshev());
  sgpp::combigrid::CombiEvaluators::Collection evaluators(
      numDimensions, sgpp::combigrid::CombiEvaluators::BSplineInterpolation(degree));
  return auxiliaryBsplineFunction(numDimensions, func, pointHierarchies, evaluators, degree);
}

std::shared_ptr<CombigridOperation> CombigridOperation::createLinearLejaBsplineInterpolation(
    size_t numDimensions, MultiFunction func, size_t degree = 3, size_t growthFactor = 2) {
  sgpp::combigrid::CombiHierarchies::Collection pointHierarchies(
      numDimensions, sgpp::combigrid::CombiHierarchies::linearLeja(growthFactor));
  sgpp::combigrid::CombiEvaluators::Collection evaluators(
      numDimensions, sgpp::combigrid::CombiEvaluators::BSplineInterpolation(degree));
  return auxiliaryBsplineFunction(numDimensions, func, pointHierarchies, evaluators, degree);
}

std::shared_ptr<CombigridOperation> CombigridOperation::createLinearL2LejaBsplineInterpolation(
    size_t numDimensions, MultiFunction func, size_t degree = 3, size_t growthFactor = 2) {
  sgpp::combigrid::CombiHierarchies::Collection pointHierarchies(
      numDimensions, sgpp::combigrid::CombiHierarchies::linearL2Leja(growthFactor));
  sgpp::combigrid::CombiEvaluators::Collection evaluators(
      numDimensions, sgpp::combigrid::CombiEvaluators::BSplineInterpolation(degree));
  return auxiliaryBsplineFunction(numDimensions, func, pointHierarchies, evaluators, degree);
}
std::shared_ptr<CombigridOperation> CombigridOperation::createExpUniformBoundaryBsplineQuadrature(
    size_t numDimensions, MultiFunction func, size_t degree = 3) {
  sgpp::combigrid::CombiHierarchies::Collection pointHierarchies(
      numDimensions, sgpp::combigrid::CombiHierarchies::expUniformBoundary());
  sgpp::combigrid::CombiEvaluators::Collection evaluators(
      numDimensions, sgpp::combigrid::CombiEvaluators::BSplineQuadrature(degree));
  return auxiliaryBsplineFunction(numDimensions, func, pointHierarchies, evaluators, degree);
}

} /* namespace combigrid */
} /* namespace sgpp*/
