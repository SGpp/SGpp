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
  numDims = weightFunctions.size();
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
  for (size_t d = 0; d < numDims; d++) {
    quadEvaluators.push_back(sgpp::combigrid::CombiEvaluators::BSplineQuadrature(
        config.degree, weightFunctions[d], numAdditionalPoints, config.bounds[2 * d],
        config.bounds[2 * d + 1], normalizeWeights));
  }
  auto quadratureOperation = std::make_shared<sgpp::combigrid::CombigridOperation>(
      pointHierarchies, quadEvaluators, levelManager, coefficientStorage, summationStrategyType);
  combigridOperation = quadratureOperation;

  scalarProducts.setWeightFunction(weightFunctions);
  scalarProducts.setBounds(config.bounds);

  currentNumGridPoints = 0;
}

// ToDo (rehmemk) tried to use addLevelsFromStructureParallel here, does not work

void BsplineStochasticCollocation::updateConfig(
    sgpp::combigrid::CombigridSurrogateModelConfiguration newConfig) {
  this->config.coefficientStorage = newConfig.coefficientStorage;
  this->config.levelStructure = newConfig.levelStructure;
  this->config.levelManager = newConfig.levelManager;

  combigridMultiOperation = createBsplineLinearCoefficientOperation(newConfig.degree, numDims,
                                                                    newConfig.coefficientStorage);
  combigridMultiOperation->getLevelManager()->addLevelsFromStructure(newConfig.levelStructure);

  size_t numAdditionalPoints = 0;
  bool normalizeWeights = false;
  sgpp::combigrid::FullGridSummationStrategyType summationStrategyType =
      sgpp::combigrid::FullGridSummationStrategyType::LINEAR;
  sgpp::combigrid::CombiEvaluators::Collection quadEvaluators(0);

  for (size_t d = 0; d < numDims; d++) {
    quadEvaluators.push_back(sgpp::combigrid::CombiEvaluators::BSplineQuadrature(
        newConfig.degree, weightFunctions[d], numAdditionalPoints, newConfig.bounds[2 * d],
        newConfig.bounds[2 * d + 1], normalizeWeights));
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

double BsplineStochasticCollocation::eval(sgpp::base::DataVector& x) {
  throw sgpp::base::application_exception("BsplineStochasticCollocation::eval - not implemented.");
}

void BsplineStochasticCollocation::eval(sgpp::base::DataMatrix& xs, sgpp::base::DataVector& res) {
  throw sgpp::base::application_exception("BsplineStochasticCollocation::eval - not implemented.");
}

double BsplineStochasticCollocation::computeMean() {
  double mean = combigridOperation->getResult();
  //  double width = 1.0;
  //  for (size_t d = 0; d < numDims; d++) {
  //    width *= (config.bounds[2 * d + 1] - config.bounds[2 * d]);
  //  }
  //  mean *= width;
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
  sgpp::base::DataVector alpha = createInterpolantOnConvertedExpUnifromBoundaryCombigird(
      grid, gridStorage, combigridMultiOperation, levelStructure);

  sgpp::base::Grid* gridptr = grid.get();
  sgpp::base::DataVector product(alpha.size());

  //  double width = 1.0;
  //  for (size_t d = 0; d < numDims; d++) {
  //    width *= (config.bounds[2 * d + 1] - config.bounds[2 * d]);
  //  }

  scalarProducts.updateGrid(gridptr);
  scalarProducts.setWeightFunction(weightFunctions);
  scalarProducts.setBounds(config.bounds);

  double variance = 0;
  if (config.degree == 1) {
    // calculate V(u) = E(u^2) - E(u)^2
    // this works for all B spline degrees
    scalarProducts.mult(alpha, product);
    double meanSquare = product.dotProduct(alpha);
    //    meanSquare *= width;
    variance = meanSquare - ev * ev;
  } else {
    // calculate V(u) = E((u-E(u))^2)
    // this is done by subtracting E(u) from the coefficient of the constant function
    // it does not work for B spline degree 1 because there is no constant function in the basis
    // (We could add a constant basis function on level 0 ()
    alpha[0] -= ev;
    scalarProducts.mult(alpha, product);
    variance = product.dotProduct(alpha);
    //    variance *= width;
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
