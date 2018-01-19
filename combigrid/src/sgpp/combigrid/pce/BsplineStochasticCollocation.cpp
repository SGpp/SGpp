// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

//#include <sgpp/combigrid/integration/GaussLegendreQuadrature.hpp>
#include <sgpp/combigrid/operation/CombigridMultiOperation.hpp>
#include <sgpp/combigrid/operation/CombigridOperation.hpp>
#include <sgpp/combigrid/operation/CombigridTensorOperation.hpp>
#include <sgpp/combigrid/pce/BsplineStochasticCollocation.hpp>
#include <sgpp/combigrid/pce/CombigridSurrogateModel.hpp>
#include <sgpp/combigrid/utils/BSplineRoutines.hpp>

#include <sgpp/base/exception/application_exception.hpp>

#include <algorithm>
#include <cmath>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace sgpp {
namespace combigrid {

// ToDo(rehmemk) implement bounds for non-unitcube domains

BsplineStochasticCollocation::BsplineStochasticCollocation(
    sgpp::combigrid::CombigridSurrogateModelConfiguration& config)
    : CombigridSurrogateModel(config),
      weightFunctions(config.weightFunctions),
      numGridPoints(0),
      computedMeanFlag(false),
      ev(0.0),
      computedVarianceFlag(false),
      var(0.0),
      coefficientStorage(config.coefficientStorage) {
  // create vector of function bases
  //  if (config.basisFunctions.size() == 0) {
  //    for (size_t idim = 0; idim < numDims; idim++) {
  //      this->config.basisFunctions.push_back(config.basisFunction);
  //    }
  //  } else if (numDims != config.basisFunctions.size()) {
  //    throw sgpp::base::application_exception(
  //        "BsplineStochasticCollocation: number of basis function do not match with the number of
  //        "
  //        "dimensions of the operation");
  //  }

  initializeOperations(config.pointHierarchies, coefficientStorage, config.levelManager);
}

BsplineStochasticCollocation::~BsplineStochasticCollocation() {}

// --------------------------------------------------------------------------------------

void BsplineStochasticCollocation::initializeOperations(
    std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
    std::shared_ptr<AbstractCombigridStorage> coefficientStorage,
    std::shared_ptr<LevelManager> levelManager) {
  // initialize interpolation operation
  sgpp::combigrid::EvaluatorConfiguration evalConfig(
      sgpp::combigrid::CombiEvaluatorTypes::Multi_BSplineInterpolation, this->config.degree);
  sgpp::combigrid::CombiEvaluators::MultiCollection evaluators(
      this->config.numDims,
      sgpp::combigrid::CombiEvaluators::createCombiMultiEvaluator(evalConfig));
  sgpp::combigrid::FullGridSummationStrategyType summationStrategyType =
      sgpp::combigrid::FullGridSummationStrategyType::LINEAR;
  auto interpolationOperation = std::make_shared<sgpp::combigrid::CombigridMultiOperation>(
      pointHierarchies, evaluators, levelManager, coefficientStorage, summationStrategyType);

  this->config.combigridMultiOperation = interpolationOperation;

  //  this->config.combigridOperation = createBsplineQuadratureCoefficientOperation(
  //      this->config.degree, this->config.numDims, levelManager, pointHierarchies,
  //      coefficientStorage);

  size_t numAdditionalPoints = 0;
  bool normalizeWeights = false;
  sgpp::combigrid::CombiEvaluators::Collection quadEvaluators(0);
  for (size_t d = 0; d < this->config.numDims; d++) {
    quadEvaluators.push_back(sgpp::combigrid::CombiEvaluators::BSplineQuadrature(
        config.degree, weightFunctions[d], numAdditionalPoints, config.bounds[2 * d],
        config.bounds[2 * d + 1], normalizeWeights));
  }
  this->config.combigridOperation = std::make_shared<sgpp::combigrid::CombigridOperation>(
      pointHierarchies, quadEvaluators, levelManager, coefficientStorage, summationStrategyType);

  numGridPoints = 0;
}

// ToDo (rehmemk) tried to use addLevelsFromStructureParallel here, does not work
void BsplineStochasticCollocation::updateConfig(
    sgpp::combigrid::CombigridSurrogateModelConfiguration newConfig) {
  this->config.coefficientStorage = newConfig.coefficientStorage;
  this->config.levelStructure = newConfig.levelStructure;

  this->config.combigridMultiOperation = createBsplineLinearCoefficientOperation(
      newConfig.degree, newConfig.numDims, newConfig.coefficientStorage);

  //  this->config.combigridOperation = createexpUniformBsplineQuadratureCoefficientOperation(
  //      newConfig.degree, newConfig.numDims, newConfig.coefficientStorage);
  size_t numAdditionalPoints = 0;
  bool normalizeWeights = false;
  sgpp::combigrid::FullGridSummationStrategyType summationStrategyType =
      sgpp::combigrid::FullGridSummationStrategyType::LINEAR;
  sgpp::combigrid::CombiEvaluators::Collection quadEvaluators(0);
  for (size_t d = 0; d < newConfig.numDims; d++) {
    quadEvaluators.push_back(sgpp::combigrid::CombiEvaluators::BSplineQuadrature(
        newConfig.degree, weightFunctions[d], numAdditionalPoints, newConfig.bounds[2 * d],
        newConfig.bounds[2 * d + 1], normalizeWeights));
  }
  this->config.combigridOperation = std::make_shared<sgpp::combigrid::CombigridOperation>(
      config.pointHierarchies, quadEvaluators, config.levelManager, newConfig.coefficientStorage,
      summationStrategyType);
  this->config.combigridOperation->getLevelManager()->addLevelsFromStructure(
      newConfig.levelStructure);
}

bool BsplineStochasticCollocation::updateStatus() {
  if (numGridPoints < config.combigridMultiOperation->numGridPoints()) {
    coefficientStorage = config.combigridMultiOperation->getStorage();
    numGridPoints = config.combigridMultiOperation->numGridPoints();
    computedMeanFlag = false;
    computedVarianceFlag = false;
    return true;
  } else {
    return false;
  }
}

double BsplineStochasticCollocation::computeMean() {
  double mean = this->config.combigridOperation->getResult();
  double width = 1.0;
  for (size_t d = 0; d < numDims; d++) {
    width *= (config.bounds[2 * d + 1] - config.bounds[2 * d]);
  }
  mean *= width;
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
  std::shared_ptr<sgpp::base::Grid> grid;
  grid.reset(sgpp::base::Grid::createNakBsplineBoundaryCombigridGrid(numDims, config.degree));
  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  auto levelStructure = this->config.levelStructure;
  convertexpUniformBoundaryCombigridToHierarchicalSparseGrid(levelStructure, gridStorage);

  // interpolate on SG
  sgpp::base::DataVector alpha = createInterpolantOnConvertedExpUnifromBoundaryCombigird(
      grid, gridStorage, this->config.combigridMultiOperation, levelStructure);

  sgpp::base::Grid* gridptr = grid.get();
  sgpp::pde::OperationMatrixLTwoDotNakBsplineBoundaryCombigrid massMatrix(gridptr, weightFunctions,
                                                                          config.bounds);
  sgpp::base::DataVector product(alpha.size(), 0);
  massMatrix.mult(alpha, product);
  double meanSquare = product.dotProduct(alpha);
  double width = 1.0;
  for (size_t d = 0; d < numDims; d++) {
    width *= (config.bounds[2 * d + 1] - config.bounds[2 * d]);
  }
  meanSquare *= width;

  if (!computedMeanFlag) {
    mean();
  }
  double variance = meanSquare - ev * ev;
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

// DEBRECATED
void BsplineStochasticCollocation::getComponentSobolIndices(
    sgpp::base::DataVector& componentSsobolIndices, bool normalized) {
  std::cout << "debrecated" << std::endl;
}
void BsplineStochasticCollocation::getTotalSobolIndices(sgpp::base::DataVector& totalSobolIndices,
                                                        bool normalized) {
  std::cout << "debrecated" << std::endl;
}

void BsplineStochasticCollocation::updateOperation(
    std::shared_ptr<sgpp::combigrid::CombigridOperation> combigridOperation) {
  std::cout << "BsplineStochasticCollocation::updateOperation debrecated" << std::endl;
}
void BsplineStochasticCollocation::updateOperation(
    std::shared_ptr<sgpp::combigrid::CombigridMultiOperation> combigridOperation) {
  std::cout << "BsplineStochasticCollocation::updateOperation (Multi) debrecated" << std::endl;
}
void BsplineStochasticCollocation::updateOperation(
    std::shared_ptr<sgpp::combigrid::CombigridTensorOperation> combigridOperation) {
  std::cout << "BsplineStochasticCollocation::updateOperation (Tensor) debrecated" << std::endl;
}

} /* namespace combigrid */
} /* namespace sgpp */
