// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/pce/PolynomialStochasticCollocation.hpp>
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
    std::shared_ptr<sgpp::combigrid::CombigridOperation> combigridOperation,
    std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D> functionBasis,
    sgpp::base::DataVector const& bounds)
    : numDims(combigridOperation->numDims()),
      combigridOperation(combigridOperation),
      combigridMultiOperation(nullptr),
      combigridTensorOperation(nullptr),
      tensorBasis(0),
      weightFunctions(0),
      bounds(bounds),
      numGridPoints(0),
      computedMeanFlag(false),
      ev(0.0),
      computedVarianceFlag(false),
      var(0.0) {
  // create vector of function bases
  for (size_t idim = 0; idim < numDims; idim++) {
    tensorBasis.push_back(functionBasis);
  }

  initializeBounds();
  initializeWeightFunctions();
  initializeTensorOperation(combigridOperation->getPointHierarchies(),
                            combigridOperation->getStorage(),
                            combigridOperation->getLevelManager());
}

PolynomialStochasticCollocation::PolynomialStochasticCollocation(
    std::shared_ptr<sgpp::combigrid::CombigridMultiOperation> combigridMultiOperation,
    std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D> functionBasis,
    sgpp::base::DataVector const& bounds)
    : numDims(combigridMultiOperation->numDims()),
      combigridOperation(nullptr),
      combigridMultiOperation(combigridMultiOperation),
      combigridTensorOperation(nullptr),
      tensorBasis(0),
      weightFunctions(0),
      bounds(bounds),
      numGridPoints(0),
      computedMeanFlag(false),
      ev(0.0),
      computedVarianceFlag(false),
      var(0.0) {
  // create vector of function bases
  for (size_t idim = 0; idim < numDims; idim++) {
    tensorBasis.push_back(functionBasis);
  }

  initializeBounds();
  initializeWeightFunctions();
  initializeTensorOperation(combigridMultiOperation->getPointHierarchies(),
                            combigridMultiOperation->getStorage(),
                            combigridMultiOperation->getLevelManager());
}

PolynomialStochasticCollocation::PolynomialStochasticCollocation(
    std::shared_ptr<sgpp::combigrid::CombigridTensorOperation> combigridTensorOperation,
    std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D> functionBasis,
    sgpp::base::DataVector const& bounds)
    : numDims(combigridTensorOperation->numDims()),
      combigridOperation(nullptr),
      combigridMultiOperation(nullptr),
      combigridTensorOperation(combigridTensorOperation),
      tensorBasis(0),
      weightFunctions(0),
      bounds(bounds),
      numGridPoints(0),
      computedMeanFlag(false),
      ev(0.0),
      computedVarianceFlag(false),
      var(0.0) {
  // create vector of function bases
  for (size_t idim = 0; idim < numDims; idim++) {
    tensorBasis.push_back(functionBasis);
  }

  initializeBounds();
  initializeWeightFunctions();
  initializeTensorOperation(combigridTensorOperation->getPointHierarchies(),
                            combigridTensorOperation->getStorage(),
                            combigridTensorOperation->getLevelManager());
}

PolynomialStochasticCollocation::PolynomialStochasticCollocation(
    std::shared_ptr<sgpp::combigrid::CombigridOperation> combigridOperation,
    sgpp::combigrid::OrthogonalBasisFunctionsCollection& tensorBasis,
    sgpp::base::DataVector const& bounds)
    : numDims(combigridOperation->numDims()),
      combigridOperation(combigridOperation),
      combigridMultiOperation(nullptr),
      combigridTensorOperation(nullptr),
      tensorBasis(tensorBasis),
      weightFunctions(0),
      bounds(bounds),
      numGridPoints(0),
      computedMeanFlag(false),
      ev(0.0),
      computedVarianceFlag(false),
      var(0.0) {
  // make sure that the number of dimensions match
  if (numDims != tensorBasis.size()) {
    throw sgpp::base::application_exception(
        "number of basis function do not match with the number of dimensions of the operation");
  }

  initializeBounds();
  initializeWeightFunctions();
  initializeTensorOperation(combigridOperation->getPointHierarchies(),
                            combigridOperation->getStorage(),
                            combigridOperation->getLevelManager());
}

PolynomialStochasticCollocation::PolynomialStochasticCollocation(
    std::shared_ptr<sgpp::combigrid::CombigridMultiOperation> combigridMultiOperation,
    sgpp::combigrid::OrthogonalBasisFunctionsCollection& tensorBasis,
    sgpp::base::DataVector const& bounds)
    : numDims(combigridMultiOperation->numDims()),
      combigridOperation(nullptr),
      combigridMultiOperation(combigridMultiOperation),
      combigridTensorOperation(nullptr),
      tensorBasis(tensorBasis),
      weightFunctions(0),
      bounds(bounds),
      numGridPoints(0),
      computedMeanFlag(false),
      ev(0.0),
      computedVarianceFlag(false),
      var(0.0) {
  // make sure that the number of dimensions match
  if (numDims != tensorBasis.size()) {
    throw sgpp::base::application_exception(
        "number of basis function do not match with the number of dimensions of the operation");
  }

  initializeBounds();
  initializeWeightFunctions();
  initializeTensorOperation(combigridMultiOperation->getPointHierarchies(),
                            combigridMultiOperation->getStorage(),
                            combigridMultiOperation->getLevelManager());
}

PolynomialStochasticCollocation::PolynomialStochasticCollocation(
    std::shared_ptr<sgpp::combigrid::CombigridTensorOperation> combigridTensorOperation,
    sgpp::combigrid::OrthogonalBasisFunctionsCollection& tensorBasis,
    sgpp::base::DataVector const& bounds)
    : numDims(combigridTensorOperation->numDims()),
      combigridOperation(nullptr),
      combigridMultiOperation(nullptr),
      combigridTensorOperation(nullptr),
      tensorBasis(tensorBasis),
      weightFunctions(0),
      bounds(bounds),
      numGridPoints(0),
      computedMeanFlag(false),
      ev(0.0),
      computedVarianceFlag(false),
      var(0.0) {
  // make sure that the number of dimensions match
  if (numDims != tensorBasis.size()) {
    throw sgpp::base::application_exception(
        "number of basis function do not match with the number of dimensions of the operation");
  }

  initializeBounds();
  initializeWeightFunctions();
  initializeTensorOperation(combigridTensorOperation->getPointHierarchies(),
                            combigridTensorOperation->getStorage(),
                            combigridTensorOperation->getLevelManager());
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

  combigridTensorOperation =
      sgpp::combigrid::CombigridTensorOperation::createOperationTensorPolynomialInterpolation(
          pointHierarchies, storage, levelManager, legendreBasis);

  numGridPoints = 0;
}

void PolynomialStochasticCollocation::initializeBounds() {
  if (bounds.size() == 0) {
    bounds.resize(2 * numDims);
    for (size_t idim = 0; idim < numDims; idim++) {
      bounds[2 * idim] = tensorBasis[idim]->lowerBound();
      bounds[2 * idim + 1] = tensorBasis[idim]->upperBound();
    }
  } else {
    if (bounds.size() != 2 * numDims) {
      throw sgpp::base::application_exception(
          "PolynomialStochasticCollocation::initializeBounds - not enough arguments for bounds "
          "specified");
    }
  }
}

void PolynomialStochasticCollocation::initializeWeightFunctions() {
  weightFunctions.clear();
  for (size_t idim = 0; idim < numDims; idim++) {
    weightFunctions.push_back(tensorBasis[idim]->getWeightFunction());
  }
}

bool PolynomialStochasticCollocation::updateStatus() {
  if (numGridPoints < combigridTensorOperation->numGridPoints()) {
    expansionCoefficients = combigridTensorOperation->getResult();
    numGridPoints = combigridTensorOperation->numGridPoints();
    computedMeanFlag = false;
    computedVarianceFlag = false;
    return true;
  } else {
    return false;
  }
}

double PolynomialStochasticCollocation::computeMean() {
  return FirstMomentNormStrategy(legendreBasis, weightFunctions, false, bounds)
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
  return VarianceNormStrategy(legendreBasis, weightFunctions, false, bounds)
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

std::shared_ptr<sgpp::combigrid::CombigridTensorOperation>
PolynomialStochasticCollocation::getCombigridTensorOperation() {
  return combigridTensorOperation;
}

void PolynomialStochasticCollocation::updateOperation(
    std::shared_ptr<sgpp::combigrid::CombigridOperation> combigridOperation) {
  initializeTensorOperation(combigridOperation->getPointHierarchies(),
                            combigridOperation->getStorage(),
                            combigridOperation->getLevelManager());
}
void PolynomialStochasticCollocation::updateOperation(
    std::shared_ptr<sgpp::combigrid::CombigridMultiOperation> combigridOperation) {
  initializeTensorOperation(combigridOperation->getPointHierarchies(),
                            combigridOperation->getStorage(),
                            combigridOperation->getLevelManager());
}
void PolynomialStochasticCollocation::updateOperation(
    std::shared_ptr<sgpp::combigrid::CombigridTensorOperation> combigridOperation) {
  initializeTensorOperation(combigridOperation->getPointHierarchies(),
                            combigridOperation->getStorage(),
                            combigridOperation->getLevelManager());
}

} /* namespace combigrid */
} /* namespace sgpp */
