// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/functions/AbstractInfiniteFunctionBasis1D.hpp>
#include <sgpp/combigrid/functions/OrthogonalPolynomialBasis1D.hpp>
#include <sgpp/combigrid/integration/GaussLegendreQuadrature.hpp>
#include <sgpp/combigrid/operation/CombigridMultiOperation.hpp>
#include <sgpp/combigrid/operation/CombigridOperation.hpp>
#include <sgpp/combigrid/operation/CombigridTensorOperation.hpp>

#include <sgpp/combigrid/operation/multidim/AveragingLevelManager.hpp>

#include <sgpp/combigrid/pce/CombigridSurrogateModel.hpp>
#include <sgpp/combigrid/pce/PolynomialStochasticCollocation.hpp>
#include <sgpp/combigrid/pce/SGppToDakota.hpp>

#include <sgpp/combigrid/algebraic/FirstMomentNormStrategy.hpp>
#include <sgpp/combigrid/algebraic/VarianceNormStrategy.hpp>

#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/tools/GaussLegendreQuadRule1D.hpp>

#include <algorithm>
#include <cmath>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace sgpp {
namespace combigrid {

PolynomialStochasticCollocation::PolynomialStochasticCollocation(
    sgpp::combigrid::CombigridSurrogateModelConfiguration& config)
    : CombigridSurrogateModel(config),
      weightFunctions(0),
      currentNumGridPoints(0),
      computedMeanFlag(false),
      ev(0.0),
      computedVarianceFlag(false),
      var(0.0) {
  // create vector of function bases
  if (config.basisFunctions.size() == 0) {
    for (size_t idim = 0; idim < numDims; idim++) {
      this->config.basisFunctions.push_back(config.basisFunction);
    }
  } else if (numDims != config.basisFunctions.size()) {
    throw sgpp::base::application_exception(
        "PolynomialStochasticCollocation: number of basis function do not match with the number of "
        "dimensions of the operation");
  }

  initializeBounds();
  initializeWeightFunctions();
  updateConfig(config);
  initializeNormStrategies();
}

PolynomialStochasticCollocation::~PolynomialStochasticCollocation() {}

// --------------------------------------------------------------------------------------

void PolynomialStochasticCollocation::initializeTensorOperation(
    std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
    std::shared_ptr<AbstractCombigridStorage> storage) {
  // create tensor operation for pce transformation
  sgpp::combigrid::OrthogonalPolynomialBasis1DConfiguration config;
  config.polyParameters.type_ = sgpp::combigrid::OrthogonalPolynomialBasisType::LEGENDRE;
  config.polyParameters.lowerBound_ = 0.0;
  config.polyParameters.upperBound_ = 1.0;
  legendreBasis = std::make_shared<sgpp::combigrid::OrthogonalPolynomialBasis1D>(config);

  this->combigridTensorOperation =
      sgpp::combigrid::CombigridTensorOperation::createOperationTensorPolynomialInterpolation(
          pointHierarchies, storage, legendreBasis);

  currentNumGridPoints = 0;
}

void PolynomialStochasticCollocation::initializeBounds() {
  if (config.bounds.size() == 0) {
    config.bounds.resize(2 * numDims);
    for (size_t idim = 0; idim < numDims; idim++) {
      config.bounds[2 * idim] = config.basisFunctions[idim]->lowerBound();
      config.bounds[2 * idim + 1] = config.basisFunctions[idim]->upperBound();
    }
  } else {
    if (config.bounds.size() != 2 * numDims) {
      throw sgpp::base::application_exception(
          "PolynomialStochasticCollocation::initializeBounds - not enough arguments for bounds "
          "specified");
    }
  }
}

void PolynomialStochasticCollocation::initializeWeightFunctions() {
  weightFunctions.clear();
  for (size_t idim = 0; idim < numDims; idim++) {
    weightFunctions.push_back(config.basisFunctions[idim]->getWeightFunction());
  }
}

void PolynomialStochasticCollocation::initializeNormStrategies() {
  firstMomentNormstrategy.reset(
      new FirstMomentNormStrategy(legendreBasis, weightFunctions, false, config.bounds));
  varianceNormStrategy.reset(
      new VarianceNormStrategy(legendreBasis, weightFunctions, false, config.bounds));
}

bool PolynomialStochasticCollocation::updateStatus() {
  if (currentNumGridPoints < combigridTensorOperation->numGridPoints()) {
    expansionCoefficients = combigridTensorOperation->getResult();
    currentNumGridPoints = combigridTensorOperation->numGridPoints();
    computedMeanFlag = false;
    computedVarianceFlag = false;
    return true;
  } else {
    return false;
  }
}

double PolynomialStochasticCollocation::computeMean() {
  return firstMomentNormstrategy->norm(expansionCoefficients);
}

double PolynomialStochasticCollocation::mean() {
  updateStatus();
  if (!computedMeanFlag) {
    ev = computeMean();
    computedMeanFlag = true;
  }
  return ev;
}

double PolynomialStochasticCollocation::computeVariance() {
  return varianceNormStrategy->norm(expansionCoefficients);
}

double PolynomialStochasticCollocation::variance() {
  updateStatus();
  if (!computedVarianceFlag) {
    var = computeVariance();
    computedVarianceFlag = true;
  }
  return var;
}

void PolynomialStochasticCollocation::getComponentSobolIndices(
    sgpp::base::DataVector& componentSsobolIndices, bool normalized) {
  throw sgpp::base::application_exception(
      "PolynomialStochasticCollocation::getComponentSobolIndices - not implemented.");
}

void PolynomialStochasticCollocation::getTotalSobolIndices(
    sgpp::base::DataVector& totalSobolIndices, bool normalized) {
  throw sgpp::base::application_exception(
      "PolynomialStochasticCollocation::getTotalSobolIndices - not implemented.");
}

void PolynomialStochasticCollocation::updateConfig(
    sgpp::combigrid::CombigridSurrogateModelConfiguration config) {
  // initialize tensor operation
  if (config.pointHierarchies.size() == numDims && config.storage != nullptr) {
    initializeTensorOperation(config.pointHierarchies, config.storage);
  }

  if (config.levelManager != nullptr) {
    combigridTensorOperation->setLevelManager(config.levelManager);
  }

  if (config.levelStructure != nullptr && combigridTensorOperation != nullptr) {
    combigridTensorOperation->getLevelManager()->addLevelsFromStructure(config.levelStructure);
  }
}

size_t PolynomialStochasticCollocation::numGridPoints() { return currentNumGridPoints; }

std::shared_ptr<LevelInfos> PolynomialStochasticCollocation::getInfoOnAddedLevels() {
  return combigridTensorOperation->getLevelManager()->getInfoOnAddedLevels();
}

} /* namespace combigrid */
} /* namespace sgpp */
