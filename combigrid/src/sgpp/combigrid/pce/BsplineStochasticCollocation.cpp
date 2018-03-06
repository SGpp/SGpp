// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/combigrid/functions/ProbabilityDensityFunction1D.hpp>
#include <sgpp/combigrid/operation/CombigridMultiOperation.hpp>
#include <sgpp/combigrid/operation/CombigridOperation.hpp>
#include <sgpp/combigrid/operation/CombigridTensorOperation.hpp>
#include <sgpp/combigrid/pce/BsplineStochasticCollocation.hpp>
#include <sgpp/combigrid/pce/CombigridSurrogateModel.hpp>
#include <sgpp/combigrid/utils/BSplineRoutines.hpp>

#include <sgpp/combigrid/utils/Stopwatch.hpp>

#include <algorithm>
#include <cmath>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace sgpp {
namespace combigrid {

BsplineStochasticCollocation::BsplineStochasticCollocation(
    sgpp::combigrid::CombigridSurrogateModelConfiguration& config)
    : CombigridSurrogateModel(config),
      weightFunctions(config.weightFunctions),
      currentNumGridPoints(0),
      computedMeanFlag(false),
      ev(0.0),
      computedVarianceFlag(false),
      var(0.0),
      coefficientStorage(config.coefficientStorage),
      scalarProducts() {
  initializeOperations(config.pointHierarchies, coefficientStorage, config.levelManager);
}

BsplineStochasticCollocation::~BsplineStochasticCollocation() {}

// combigridMultiOperation is a BSplineInterpolation
// combigridOperation is a BSplineQuadrature
void BsplineStochasticCollocation::initializeOperations(
    std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
    std::shared_ptr<AbstractCombigridStorage> coefficientStorage,
    std::shared_ptr<LevelManager> levelManager) {
  customWeightFunction = true;
  numDims = weightFunctions.size();
  if (numDims == 0) {
    numDims = config.numDimensions;
    customWeightFunction = false;
  }
  // initialize interpolation operation
  sgpp::combigrid::EvaluatorConfiguration evalConfig(
      sgpp::combigrid::CombiEvaluatorTypes::Multi_BSplineInterpolation, this->config.degree);
  sgpp::combigrid::CombiEvaluators::MultiCollection evaluators(
      numDims, sgpp::combigrid::CombiEvaluators::createCombiMultiEvaluator(evalConfig));
  sgpp::combigrid::FullGridSummationStrategyType summationStrategyType =
      sgpp::combigrid::FullGridSummationStrategyType::LINEAR;
  auto interpolationOperation = std::make_shared<sgpp::combigrid::CombigridMultiOperation>(
      pointHierarchies, evaluators, levelManager, coefficientStorage, summationStrategyType);

  combigridMultiOperation = interpolationOperation;

  size_t numAdditionalPoints = 0;
  bool normalizeWeights = false;
  sgpp::combigrid::CombiEvaluators::Collection quadEvaluators(0);
  if (customWeightFunction) {
    for (size_t d = 0; d < numDims; d++) {
      quadEvaluators.push_back(sgpp::combigrid::CombiEvaluators::BSplineQuadrature(
          config.degree, weightFunctions[d], numAdditionalPoints, config.bounds[2 * d],
          config.bounds[2 * d + 1], normalizeWeights));
    }

    scalarProducts.setWeightFunction(weightFunctions);
    scalarProducts.setBounds(config.bounds);
  } else {
    for (size_t d = 0; d < numDims; d++) {
      quadEvaluators.push_back(sgpp::combigrid::CombiEvaluators::BSplineQuadrature(config.degree));
    }
  }
  auto quadratureOperation = std::make_shared<sgpp::combigrid::CombigridOperation>(
      pointHierarchies, quadEvaluators, levelManager, coefficientStorage, summationStrategyType);
  combigridOperation = quadratureOperation;

  currentNumGridPoints = 0;
}

// ToDo (rehmemk) tried to use addLevelsFromStructureParallel here, does not work

void BsplineStochasticCollocation::updateConfig(
    sgpp::combigrid::CombigridSurrogateModelConfiguration newConfig) {
  this->config.coefficientStorage = newConfig.coefficientStorage;
  this->config.levelStructure = newConfig.levelStructure;

  combigridMultiOperation = CombigridMultiOperation::createBsplineLinearCoefficientOperation(
      newConfig.degree, numDims, newConfig.coefficientStorage);
  combigridMultiOperation->getLevelManager()->addLevelsFromStructure(newConfig.levelStructure);

  size_t numAdditionalPoints = 0;
  bool normalizeWeights = false;
  sgpp::combigrid::FullGridSummationStrategyType summationStrategyType =
      sgpp::combigrid::FullGridSummationStrategyType::LINEAR;
  sgpp::combigrid::CombiEvaluators::Collection quadEvaluators(0);

  if (customWeightFunction) {
    for (size_t d = 0; d < numDims; d++) {
      quadEvaluators.push_back(sgpp::combigrid::CombiEvaluators::BSplineQuadrature(
          config.degree, weightFunctions[d], numAdditionalPoints, config.bounds[2 * d],
          config.bounds[2 * d + 1], normalizeWeights));
    }
  } else {
    for (size_t d = 0; d < numDims; d++) {
      quadEvaluators.push_back(sgpp::combigrid::CombiEvaluators::BSplineQuadrature(config.degree));
    }
  }
  combigridOperation = std::make_shared<sgpp::combigrid::CombigridOperation>(
      config.pointHierarchies, quadEvaluators, config.levelManager, newConfig.coefficientStorage,
      summationStrategyType);
  combigridOperation->getLevelManager()->addLevelsFromStructure(newConfig.levelStructure);

  computedMeanFlag = false;
  computedVarianceFlag = false;
}

bool BsplineStochasticCollocation::updateStatus() {
  if (currentNumGridPoints < combigridMultiOperation->numGridPoints()) {
    coefficientStorage = combigridMultiOperation->getStorage();
    currentNumGridPoints = combigridMultiOperation->numGridPoints();
    computedMeanFlag = false;
    computedVarianceFlag = false;
    return true;
  } else {
    return false;
  }
}

void BsplineStochasticCollocation::eval(sgpp::base::DataMatrix& xs, sgpp::base::DataVector& res) {
  combigridMultiOperation->setParameters(xs);
  combigridMultiOperation->getLevelManager()->addLevelsFromStructure(config.levelStructure);
  res = combigridMultiOperation->getResult();
}

double BsplineStochasticCollocation::eval(sgpp::base::DataVector& x) {
  sgpp::base::DataMatrix xs(x.size(), 0);
  xs.appendCol(x);
  sgpp::base::DataVector res;
  this->eval(xs, res);
  return res[0];
}

double BsplineStochasticCollocation::computeMean() {
  double mean = combigridOperation->getResult();
  return mean;
}

double BsplineStochasticCollocation::mean() {
  updateStatus();
  if (!computedMeanFlag) {
    ev = computeMean();
    computedMeanFlag = true;
  }
  return ev;
}

double BsplineStochasticCollocation::computeVariance() {
  if (!computedMeanFlag) {
    mean();
  }

  std::shared_ptr<sgpp::base::Grid> grid;
  grid.reset(sgpp::base::Grid::createNakBsplineBoundaryCombigridGrid(numDims, config.degree));
  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  auto levelStructure = this->config.levelStructure;
  convertexpUniformBoundaryCombigridToHierarchicalSparseGrid(levelStructure, gridStorage);

  // interpolate on SG
  sgpp::base::DataVector alpha =
      calculateInterpolationCoefficientsForConvertedExpUniformBoundaryCombigird(
          grid, gridStorage, combigridMultiOperation, levelStructure);

  sgpp::base::Grid* gridptr = grid.get();
  sgpp::base::DataVector product(alpha.size());

  scalarProducts.updateGrid(gridptr);
  // scalarProducts.setWeightFunction(weightFunctions);
  // scalarProducts.setBounds(config.bounds);

  double variance = 0;
  if (config.degree == 1) {
    // calculate V(u) = E(u^2) - E(u)^2
    // this works for all B spline degrees but may be instable
    scalarProducts.mult(alpha, product);
    double meanSquare = product.dotProduct(alpha);
    variance = meanSquare - ev * ev;
  } else {
    // calculate V(u) = E((u-E(u))^2)
    // this is done by subtracting E(u) from the coefficient of the constant function
    // it does not work for B spline degree 1 because there is no constant function in the basis
    // (We could add a constant basis function on level 0)
    alpha[0] -= ev;
    scalarProducts.mult(alpha, product);

    variance = product.dotProduct(alpha);
  }

  return variance;
}

double BsplineStochasticCollocation::variance() {
  updateStatus();
  if (!computedVarianceFlag) {
    var = computeVariance();
    computedVarianceFlag = true;
  }
  return var;
}

void BsplineStochasticCollocation::getComponentSobolIndices(
    sgpp::base::DataVector& componentSsobolIndices, bool normalized) {
  throw sgpp::base::application_exception(
      "BsplineStochasticCollocation::getComponentSobolIndices - not implemented.");
}
void BsplineStochasticCollocation::getTotalSobolIndices(sgpp::base::DataVector& totalSobolIndices,
                                                        bool normalized) {
  throw sgpp::base::application_exception(
      "PolynomialStochasticCollocation::getTotalSobolIndices - not implemented.");
}

size_t BsplineStochasticCollocation::numGridPoints() { return currentNumGridPoints; }

std::shared_ptr<LevelInfos> BsplineStochasticCollocation::getInfoOnAddedLevels() {
  return combigridOperation->getLevelManager()->getInfoOnAddedLevels();
}

} /* namespace combigrid */
} /* namespace sgpp */
