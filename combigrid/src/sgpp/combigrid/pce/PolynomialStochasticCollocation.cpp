// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/pce/PolynomialStochasticCollocation.hpp>
#include <sgpp/combigrid/pce/CombigridSurrogateModel.hpp>
#include <sgpp/combigrid/pce/SGppToDakota.hpp>
#include <sgpp/combigrid/operation/CombigridOperation.hpp>
#include <sgpp/combigrid/operation/CombigridMultiOperation.hpp>
#include <sgpp/combigrid/operation/CombigridTensorOperation.hpp>
#include <sgpp/combigrid/functions/AbstractInfiniteFunctionBasis1D.hpp>
#include <sgpp/combigrid/functions/OrthogonalPolynomialBasis1D.hpp>
#include <sgpp/combigrid/integration/GaussLegendreQuadrature.hpp>

#include <sgpp/combigrid/algebraic/FirstMomentNormStrategy.hpp>
#include <sgpp/combigrid/algebraic/VarianceNormStrategy.hpp>

#include <sgpp/base/tools/GaussLegendreQuadRule1D.hpp>
#include <sgpp/base/exception/application_exception.hpp>

#include <algorithm>
#include <vector>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace sgpp {
namespace combigrid {

PolynomialStochasticCollocation::PolynomialStochasticCollocation(
    sgpp::combigrid::CombigridSurrogateModelConfiguration& config)
    : CombigridSurrogateModel(config),
      weightFunctions(0),
      numGridPoints(0),
      computedMeanFlag(false),
      ev(0.0),
      computedVarianceFlag(false),
      var(0.0) {
  // create vector of function bases
  if (config.basisFunctions.size() == 0) {
    for (size_t idim = 0; idim < numDims; idim++) {
      config.basisFunctions.push_back(config.basisFunction);
    }
  } else if (numDims != config.basisFunctions.size()) {
    throw sgpp::base::application_exception(
        "PolynomialStochasticCollocation: number of basis function do not match with the number of "
        "dimensions of the operation");
  }

  initializeBounds();
  initializeWeightFunctions();

  if (config.combigridOperation != nullptr) {
    initializeTensorOperation(config.combigridOperation->getPointHierarchies(),
                              config.combigridOperation->getStorage(),
                              config.combigridOperation->getLevelManager());
  } else if (config.combigridMultiOperation != nullptr) {
    initializeTensorOperation(config.combigridMultiOperation->getPointHierarchies(),
                              config.combigridMultiOperation->getStorage(),
                              config.combigridMultiOperation->getLevelManager());
  } else if (config.combigridTensorOperation != nullptr) {
    initializeTensorOperation(config.combigridTensorOperation->getPointHierarchies(),
                              config.combigridTensorOperation->getStorage(),
                              config.combigridTensorOperation->getLevelManager());
  } else {
    throw sgpp::base::application_exception(
        "PolynomialStochasticCollocation: no operation is set in surrogate model config");
  }
}

PolynomialStochasticCollocation::~PolynomialStochasticCollocation() {}

// --------------------------------------------------------------------------------------

void PolynomialStochasticCollocation::initializeTensorOperation(
    std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
    std::shared_ptr<AbstractCombigridStorage> storage, std::shared_ptr<LevelManager> levelManager) {
  // create tensor operation for pce transformation
  sgpp::combigrid::OrthogonalPolynomialBasis1DConfiguration config;
  config.polyParameters.type_ = sgpp::combigrid::OrthogonalPolynomialBasisType::LEGENDRE;
  config.polyParameters.lowerBound_ = 0.0;
  config.polyParameters.upperBound_ = 1.0;
  legendreBasis = std::make_shared<sgpp::combigrid::OrthogonalPolynomialBasis1D>(config);

  this->config.combigridTensorOperation =
      sgpp::combigrid::CombigridTensorOperation::createOperationTensorPolynomialInterpolation(
          pointHierarchies, storage, levelManager, legendreBasis);

  numGridPoints = 0;
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

bool PolynomialStochasticCollocation::updateStatus() {
  if (numGridPoints < config.combigridTensorOperation->numGridPoints()) {
    expansionCoefficients = config.combigridTensorOperation->getResult();
    numGridPoints = config.combigridTensorOperation->numGridPoints();
    computedMeanFlag = false;
    computedVarianceFlag = false;
    return true;
  } else {
    return false;
  }
}

double PolynomialStochasticCollocation::computeMean() {
  return FirstMomentNormStrategy(legendreBasis, weightFunctions, false, config.bounds)
      .norm(expansionCoefficients);
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
  return VarianceNormStrategy(legendreBasis, weightFunctions, false, config.bounds)
      .norm(expansionCoefficients);
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

void PolynomialStochasticCollocation::updateOperation(
    std::shared_ptr<sgpp::combigrid::CombigridOperation> combigridOperation) {
  if (combigridOperation != nullptr) {
    this->config.combigridOperation = combigridOperation;
    initializeTensorOperation(combigridOperation->getPointHierarchies(),
                              combigridOperation->getStorage(),
                              combigridOperation->getLevelManager());
  }
}
void PolynomialStochasticCollocation::updateOperation(
    std::shared_ptr<sgpp::combigrid::CombigridMultiOperation> combigridOperation) {
  if (combigridOperation != nullptr) {
    this->config.combigridMultiOperation = combigridOperation;
    initializeTensorOperation(combigridOperation->getPointHierarchies(),
                              combigridOperation->getStorage(),
                              combigridOperation->getLevelManager());
  }
}
void PolynomialStochasticCollocation::updateOperation(
    std::shared_ptr<sgpp::combigrid::CombigridTensorOperation> combigridOperation) {
  if (combigridOperation != nullptr) {
    this->config.combigridTensorOperation = combigridOperation;
    initializeTensorOperation(combigridOperation->getPointHierarchies(),
                              combigridOperation->getStorage(),
                              combigridOperation->getLevelManager());
  }
}

} /* namespace combigrid */
} /* namespace sgpp */
