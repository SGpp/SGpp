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
      functionBases(0),
      bounds(bounds),
      numGridPoints(0),
      computedMeanFlag(false),
      ev(0.0),
      computedVarianceFlag(false),
      var(0.0) {
  // create vector of function bases
  for (size_t idim = 0; idim < numDims; idim++) {
    functionBases.push_back(functionBasis);
  }

  initializeBounds();
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
      functionBases(0),
      bounds(bounds),
      numGridPoints(0),
      computedMeanFlag(false),
      ev(0.0),
      computedVarianceFlag(false),
      var(0.0) {
  // create vector of function bases
  for (size_t idim = 0; idim < numDims; idim++) {
    functionBases.push_back(functionBasis);
  }

  initializeBounds();
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
      functionBases(0),
      bounds(bounds),
      numGridPoints(0),
      computedMeanFlag(false),
      ev(0.0),
      computedVarianceFlag(false),
      var(0.0) {
  // create vector of function bases
  for (size_t idim = 0; idim < numDims; idim++) {
    functionBases.push_back(functionBasis);
  }

  initializeBounds();
  initializeTensorOperation(combigridTensorOperation->getPointHierarchies(),
                            combigridTensorOperation->getStorage(),
                            combigridTensorOperation->getLevelManager());
}

PolynomialStochasticCollocation::PolynomialStochasticCollocation(
    std::shared_ptr<sgpp::combigrid::CombigridOperation> combigridOperation,
    std::vector<std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D>>& functionBases,
    sgpp::base::DataVector const& bounds)
    : numDims(combigridOperation->numDims()),
      combigridOperation(combigridOperation),
      combigridMultiOperation(nullptr),
      combigridTensorOperation(nullptr),
      functionBases(functionBases),
      bounds(bounds),
      numGridPoints(0),
      computedMeanFlag(false),
      ev(0.0),
      computedVarianceFlag(false),
      var(0.0) {
  // make sure that the number of dimensions match
  if (numDims != functionBases.size()) {
    throw sgpp::base::application_exception(
        "number of basis function do not match with the number of dimensions of the operation");
  }

  initializeBounds();
  initializeTensorOperation(combigridOperation->getPointHierarchies(),
                            combigridOperation->getStorage(),
                            combigridOperation->getLevelManager());
}

PolynomialStochasticCollocation::PolynomialStochasticCollocation(
    std::shared_ptr<sgpp::combigrid::CombigridMultiOperation> combigridMultiOperation,
    std::vector<std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D>>& functionBases,
    sgpp::base::DataVector const& bounds)
    : numDims(combigridMultiOperation->numDims()),
      combigridOperation(nullptr),
      combigridMultiOperation(combigridMultiOperation),
      combigridTensorOperation(nullptr),
      functionBases(functionBases),
      bounds(bounds),
      numGridPoints(0),
      computedMeanFlag(false),
      ev(0.0),
      computedVarianceFlag(false),
      var(0.0) {
  // make sure that the number of dimensions match
  if (numDims != functionBases.size()) {
    throw sgpp::base::application_exception(
        "number of basis function do not match with the number of dimensions of the operation");
  }

  initializeBounds();
  initializeTensorOperation(combigridMultiOperation->getPointHierarchies(),
                            combigridMultiOperation->getStorage(),
                            combigridMultiOperation->getLevelManager());
}

PolynomialStochasticCollocation::PolynomialStochasticCollocation(
    std::shared_ptr<sgpp::combigrid::CombigridTensorOperation> combigridTensorOperation,
    std::vector<std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D>>& functionBases,
    sgpp::base::DataVector const& bounds)
    : numDims(combigridTensorOperation->numDims()),
      combigridOperation(nullptr),
      combigridMultiOperation(nullptr),
      combigridTensorOperation(nullptr),
      functionBases(functionBases),
      bounds(bounds),
      numGridPoints(0),
      computedMeanFlag(false),
      ev(0.0),
      computedVarianceFlag(false),
      var(0.0) {
  // make sure that the number of dimensions match
  if (numDims != functionBases.size()) {
    throw sgpp::base::application_exception(
        "number of basis function do not match with the number of dimensions of the operation");
  }

  initializeBounds();
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
      bounds[2 * idim] = functionBases[idim]->lowerBound();
      bounds[2 * idim + 1] = functionBases[idim]->upperBound();
    }
  } else {
    if (bounds.size() != 2 * numDims) {
      throw sgpp::base::application_exception(
          "PolynomialStochasticCollocation::initializeBounds - not enough arguments for bounds "
          "specified");
    }
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

double PolynomialStochasticCollocation::quad(MultiIndex i, MultiIndex j) {
  double ans = 1.0;

  // performing Gauss-Legendre integration
  GaussLegendreQuadrature gaussLegendreQuadrature;

  // Gauss quadrature in each dimension
  for (size_t idim = 0; idim < i.size(); idim++) {
    size_t degree_i = i[idim];
    size_t degree_j = j[idim];
    auto functionBasis = functionBases[idim];
    size_t incrementQuadraturePoints = functionBasis->numAdditionalQuadraturePoints();
    size_t numGaussPoints = (degree_i + degree_j + 3) / 2;

    auto func = [&functionBasis, &degree_i, &degree_j, &idim, this](double x_unit, double x_prob) {
      return this->legendreBasis->evaluate(degree_i, x_unit) *
             this->legendreBasis->evaluate(degree_j, x_unit) * functionBasis->pdf(x_prob);
    };

    double a = bounds[2 * idim], b = bounds[2 * idim + 1];
    ans *= GaussLegendreQuadrature::evaluate_iteratively(func, a, b, numGaussPoints,
                                                         incrementQuadraturePoints);
  }
  return ans;
}

double PolynomialStochasticCollocation::quad(MultiIndex i) {
  double ans = 1.0;

  // performing Gauss-Legendre integration
  GaussLegendreQuadrature gaussLegendreQuadrature;

  // Gauss quadrature in each dimension
  for (size_t idim = 0; idim < i.size(); idim++) {
    size_t degree_i = i[idim];
    auto functionBasis = functionBases[idim];
    size_t incrementQuadraturePoints = functionBasis->numAdditionalQuadraturePoints();
    size_t numGaussPoints = (degree_i + 2) / 2;

    auto func = [&functionBasis, &degree_i, &idim, this](double x_unit, double x_prob) {
      return this->legendreBasis->evaluate(degree_i, x_unit) * functionBasis->pdf(x_prob);
    };

    double a = bounds[2 * idim], b = bounds[2 * idim + 1];
    ans *= GaussLegendreQuadrature::evaluate_iteratively(func, a, b, numGaussPoints,
                                                         incrementQuadraturePoints);
  }
  return ans;
}

double PolynomialStochasticCollocation::computeMean() {
  // compute mass matrix and corresponding coefficient vector
  auto multiIndices_i = expansionCoefficients.getValues();
  auto it_i = multiIndices_i->getStoredDataIterator();
  sgpp::base::DataVector m(numGridPoints);
  sgpp::base::DataVector coeffs(numGridPoints);

  size_t i = 0;
  while (it_i->isValid() && i < numGridPoints) {
    auto ix = it_i->getMultiIndex();

    double value = 0.0;
    auto it_value = innerProducts.find(ix);
    if (it_value != innerProducts.end()) {
      value = it_value->second;
    } else {
      value = quad(ix);
      innerProducts[ix] = value;
    }
    m[i] = value;

    coeffs[i] = it_i->value().value();
    it_i->moveToNext();
    i += 1;
  }

  // compute mean: m^T coeffs
  double ans = 0.0;
  for (size_t i = 0; i < m.getSize(); i++) {
    ans += m[i] * coeffs[i];
  }

  return ans;
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
  // prepare the mean
  double ev = mean();

  // compute mass matrix and corresponding coefficient vector
  auto multiIndices_i = expansionCoefficients.getValues();
  auto multiIndices_j = expansionCoefficients.getValues();
  sgpp::base::DataMatrix M(numGridPoints, numGridPoints);
  sgpp::base::DataVector coeffs(numGridPoints);

  auto it_i = multiIndices_i->getStoredDataIterator();
  size_t i = 0;
  while (it_i->isValid() && i < numGridPoints) {
    MultiIndex ix = it_i->getMultiIndex();

    auto it_j = multiIndices_j->getStoredDataIterator();
    size_t j = 0;
    while (it_j->isValid() && j < numGridPoints) {
      if (j >= i) {
        MultiIndex jx = it_j->getMultiIndex();

        // compute the inner product and store it
        MultiIndex kx;
        joinMultiIndices(ix, jx, kx);
        double innerProduct = 0.0;
        auto it_value = innerProducts.find(kx);
        if (it_value != innerProducts.end()) {
          innerProduct = it_value->second;
        } else {
          innerProduct = quad(ix, jx);
          innerProducts[kx] = innerProduct;
        }

        M.set(i, j, innerProduct);
        M.set(j, i, innerProduct);

        // store the coefficient in the correct order
        if (i == 0) {
          coeffs[j] = it_j->value().value();
          if (jx == MultiIndex(numDims, 0)) {
            coeffs[j] -= ev;
          }
        }
      }

      it_j->moveToNext();
      j += 1;
    }

    it_i->moveToNext();
    i += 1;
  }

  // compute variance: coeffs^T M coeffs
  sgpp::base::DataVector result(coeffs.getSize());
  M.mult(coeffs, result);

  double ans = 0.0;
  for (size_t i = 0; i < M.getNrows(); i++) {
    ans += coeffs[i] * result[i];
  }

  return ans;
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

void PolynomialStochasticCollocation::joinMultiIndices(MultiIndex& ix, MultiIndex& jx,
                                                       MultiIndex& kx) {
  kx.insert(kx.end(), ix.begin(), ix.end());
  kx.insert(kx.end(), jx.begin(), jx.end());
}

} /* namespace combigrid */
} /* namespace sgpp */
