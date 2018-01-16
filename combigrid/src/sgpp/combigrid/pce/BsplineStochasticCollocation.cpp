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

BsplineStochasticCollocation::BsplineStochasticCollocation(
    sgpp::combigrid::CombigridSurrogateModelConfiguration& config)
    : CombigridSurrogateModel(config),
      weightFunctions(0),
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

  initializeBounds();
  initializeWeightFunctions();
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

  this->config.combigridOperation = createBsplineQuadratureCoefficientOperation(
      this->config.degree, this->config.numDims, levelManager, pointHierarchies,
      coefficientStorage);

  numGridPoints = 0;
}

// ToDo (rehmemk) tried to use addLEvelsFromStructureParallel here, does not work (segmentation
// fault (core dumped))
void BsplineStochasticCollocation::updateConfig(
    sgpp::combigrid::CombigridSurrogateModelConfiguration newConfig) {
  this->config.coefficientStorage = newConfig.coefficientStorage;
  this->config.levelManager = newConfig.levelManager;
  initializeOperations(newConfig.pointHierarchies, newConfig.coefficientStorage, newConfig.levelManager);

  std::cout << "initialized new operations" << std::endl;

  std::shared_ptr<sgpp::combigrid::TreeStorage<uint8_t>> levelStructure =
      newConfig.levelManager->getLevelStructure();

  this->config.combigridMultiOperation->getLevelManager()->addLevelsFromStructure(levelStructure);
  this->config.combigridOperation->getLevelManager()->addLevelsFromStructure(levelStructure);

  std::cout << "added levels from structure" << std::endl;

  std::vector<double> res =
      calculateBsplineMeanAndVariance(newConfig.levelManager->getLevelStructure(), newConfig.numDims,
                                      newConfig.degree, newConfig.coefficientStorage);
  std::cout << "var = " << res[1] << std::endl;
}

void BsplineStochasticCollocation::initializeWeightFunctions() {
  //  weightFunctions.clear();
  //  for (size_t idim = 0; idim < numDims; idim++) {
  //    weightFunctions.push_back(config.basisFunctions[idim]->getWeightFunction());
  //  }
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
  auto levelStructure = this->config.levelManager->getLevelStructure();
  std::cout << "a" << std::endl;
  std::vector<double> res = calculateBsplineMeanAndVariance(
      levelStructure, config.numDims, config.degree, config.coefficientStorage);
  std::cout << "b" << std::endl;
  return res[1];

  //  std::shared_ptr<sgpp::base::Grid> grid;
  //  grid.reset(sgpp::base::Grid::createNakBsplineBoundaryCombigridGrid(numDims, config.degree));
  //  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  //  auto levelStructure =
  //      this->config.combigridMultiOperation->getLevelManager()->getLevelStructure();
  //  convertexpUniformBoundaryCombigridToHierarchicalSparseGrid(levelStructure, gridStorage);
  //
  //  // interpolate on SG
  //  sgpp::base::DataMatrix interpolParams(numDims, gridStorage.getSize());
  //  for (size_t i = 0; i < gridStorage.getSize(); i++) {
  //    sgpp::base::GridPoint& gp = gridStorage.getPoint(i);
  //    sgpp::base::DataVector p(gridStorage.getDimension(), 0.0);
  //    for (size_t j = 0; j < gridStorage.getDimension(); j++) {
  //      p[j] = gp.getStandardCoordinate(j);
  //    }
  //    interpolParams.setColumn(i, p);
  //  }
  //
  //  // obtain function values from combigrid surrogate
  //  this->config.combigridMultiOperation->setParameters(interpolParams);
  //  //
  //  this->config.combigridMultiOperation->getLevelManager()->addLevelsFromStructure(levelStructure);
  //  sgpp::base::DataVector f_values = this->config.combigridMultiOperation->getResult();
  //
  //  sgpp::optimization::Printer::getInstance().setVerbosity(-1);
  //  sgpp::optimization::HierarchisationSLE hierSLE(*grid);
  //  sgpp::optimization::sle_solver::UMFPACK sleSolver;
  //  sgpp::base::DataVector alpha(grid->getSize());
  //  std::cout << "a" << std::endl;
  //  if (!sleSolver.solve(hierSLE, f_values, alpha)) {
  //    std::cout << "Solving failed!" << std::endl;
  //  }
  //  std::cout << "b" << std::endl;
  //
  //  sgpp::base::Grid* gridptr = grid.get();
  //  sgpp::pde::OperationMatrixLTwoDotNakBsplineBoundaryCombigrid massMatrix(gridptr);
  //  sgpp::base::DataVector product(alpha.size(), 0);
  //  massMatrix.mult(alpha, product);
  //  double meanSquare = product.dotProduct(alpha);
  //  if (!computedMeanFlag) {
  //    computeMean();
  //  }
  //  double variance = meanSquare - ev * ev;
  //  return variance;
}

double BsplineStochasticCollocation::variance() {
  updateStatus();
  if (!computedVarianceFlag) {
    var = computeVariance();
    computedVarianceFlag = true;
  }
  return var;
}

void BsplineStochasticCollocation::initializeBounds() {
  // implement this!
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
