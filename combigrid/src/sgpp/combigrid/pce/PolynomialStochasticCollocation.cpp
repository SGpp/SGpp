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
#include <sgpp/base/tools/GaussLegendreQuadRule1D.hpp>
#include <sgpp/base/exception/application_exception.hpp>

#include <algorithm>
#include <vector>

namespace sgpp {
namespace combigrid {

PolynomialStochasticCollocation::PolynomialStochasticCollocation(
    std::shared_ptr<sgpp::combigrid::CombigridOperation> combigridOperation,
    std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D> functionBasis,
    sgpp::base::DataVector& bounds)
    : numDims(combigridOperation->numDims()),
      combigridOperation(combigridOperation),
      combigridMultiOperation(nullptr),
      combigridTensorOperation(nullptr),
      functionBases(0),
      trans(bounds),
      numGridPoints(0) {
  // create vector of function bases
  for (size_t idim = 0; idim < numDims; idim++) {
    functionBases.push_back(functionBasis);
  }

  initializeTensorOperation(combigridOperation->getPointHierarchies(),
                            combigridOperation->getStorage(),
                            combigridOperation->getLevelManager());
}

PolynomialStochasticCollocation::PolynomialStochasticCollocation(
    std::shared_ptr<sgpp::combigrid::CombigridMultiOperation> combigridMultiOperation,
    std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D> functionBasis,
    sgpp::base::DataVector& bounds)
    : numDims(combigridMultiOperation->numDims()),
      combigridOperation(nullptr),
      combigridMultiOperation(combigridMultiOperation),
      combigridTensorOperation(nullptr),
      functionBases(0),
      trans(bounds),
      numGridPoints(0) {
  // create vector of function bases
  for (size_t idim = 0; idim < numDims; idim++) {
    functionBases.push_back(functionBasis);
  }

  // create tensor operation for basis transformation
  initializeTensorOperation(combigridMultiOperation->getPointHierarchies(),
                            combigridMultiOperation->getStorage(),
                            combigridMultiOperation->getLevelManager());
}

PolynomialStochasticCollocation::PolynomialStochasticCollocation(
    std::shared_ptr<sgpp::combigrid::CombigridTensorOperation> combigridTensorOperation,
    std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D> functionBasis,
    sgpp::base::DataVector& bounds)
    : numDims(combigridTensorOperation->numDims()),
      combigridOperation(nullptr),
      combigridMultiOperation(nullptr),
      combigridTensorOperation(combigridTensorOperation),
      functionBases(0),
      trans(bounds),
      numGridPoints(0) {
  // create vector of function bases
  for (size_t idim = 0; idim < numDims; idim++) {
    functionBases.push_back(functionBasis);
  }

  initializeTensorOperation(combigridTensorOperation->getPointHierarchies(),
                            combigridTensorOperation->getStorage(),
                            combigridTensorOperation->getLevelManager());
}

PolynomialStochasticCollocation::PolynomialStochasticCollocation(
    std::shared_ptr<sgpp::combigrid::CombigridOperation> combigridOperation,
    std::vector<std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D>>& functionBases,
    sgpp::base::DataVector& bounds)
    : numDims(combigridMultiOperation->numDims()),
      combigridOperation(combigridOperation),
      combigridMultiOperation(nullptr),
      combigridTensorOperation(nullptr),
      functionBases(functionBases),
      trans(bounds),
      numGridPoints(0) {
  // make sure that the number of dimensions match
  if (numDims != functionBases.size()) {
    throw sgpp::base::application_exception(
        "number of basis function do not match with the number of dimensions of the operation");
  }

  // create tensor operation for basis transformation
  initializeTensorOperation(combigridOperation->getPointHierarchies(),
                            combigridOperation->getStorage(),
                            combigridOperation->getLevelManager());
}

PolynomialStochasticCollocation::PolynomialStochasticCollocation(
    std::shared_ptr<sgpp::combigrid::CombigridMultiOperation> combigridMultiOperation,
    std::vector<std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D>>& functionBases,
    sgpp::base::DataVector& bounds)
    : numDims(combigridMultiOperation->numDims()),
      combigridOperation(nullptr),
      combigridMultiOperation(combigridMultiOperation),
      combigridTensorOperation(nullptr),
      functionBases(functionBases),
      trans(bounds),
      numGridPoints(0) {
  // make sure that the number of dimensions match
  if (numDims != functionBases.size()) {
    throw sgpp::base::application_exception(
        "number of basis function do not match with the number of dimensions of the operation");
  }

  initializeTensorOperation(combigridMultiOperation->getPointHierarchies(),
                            combigridMultiOperation->getStorage(),
                            combigridMultiOperation->getLevelManager());
}

PolynomialStochasticCollocation::PolynomialStochasticCollocation(
    std::shared_ptr<sgpp::combigrid::CombigridTensorOperation> combigridTensorOperation,
    std::vector<std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D>>& functionBases,
    sgpp::base::DataVector& bounds)
    : numDims(combigridTensorOperation->numDims()),
      combigridOperation(nullptr),
      combigridMultiOperation(nullptr),
      combigridTensorOperation(nullptr),
      functionBases(functionBases),
      trans(bounds),
      numGridPoints(0) {
  // make sure that the number of dimensions match
  if (numDims != functionBases.size()) {
    throw sgpp::base::application_exception(
        "number of basis function do not match with the number of dimensions of the operation");
  }

  initializeTensorOperation(combigridTensorOperation->getPointHierarchies(),
                            combigridTensorOperation->getStorage(),
                            combigridTensorOperation->getLevelManager());
}

#ifdef USE_DAKOTA
PolynomialStochasticCollocation::PolynomialStochasticCollocation(
    std::shared_ptr<sgpp::combigrid::CombigridOperation> combigridOperation,
    std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D> functionBasis)
    : numDims(combigridOperation->numDims()),
      combigridOperation(combigridOperation),
      combigridMultiOperation(nullptr),
      combigridTensorOperation(nullptr),
      functionBases(0),
      numGridPoints(0) {
  // create vector of function bases
  for (size_t idim = 0; idim < numDims; idim++) {
    functionBases.push_back(functionBasis);
  }

  initializeLinearTransformation();
  initializeTensorOperation(combigridOperation->getPointHierarchies(),
                            combigridOperation->getStorage(),
                            combigridOperation->getLevelManager());
}

PolynomialStochasticCollocation::PolynomialStochasticCollocation(
    std::shared_ptr<sgpp::combigrid::CombigridMultiOperation> combigridMultiOperation,
    std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D> functionBasis)
    : numDims(combigridMultiOperation->numDims()),
      combigridOperation(nullptr),
      combigridMultiOperation(combigridMultiOperation),
      combigridTensorOperation(nullptr),
      functionBases(0),
      numGridPoints(0) {
  // create vector of function bases
  for (size_t idim = 0; idim < numDims; idim++) {
    functionBases.push_back(functionBasis);
  }

  initializeLinearTransformation();
  initializeTensorOperation(combigridMultiOperation->getPointHierarchies(),
                            combigridMultiOperation->getStorage(),
                            combigridMultiOperation->getLevelManager());
}

PolynomialStochasticCollocation::PolynomialStochasticCollocation(
    std::shared_ptr<sgpp::combigrid::CombigridTensorOperation> combigridTensorOperation,
    std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D> functionBasis)
    : numDims(combigridTensorOperation->numDims()),
      combigridOperation(nullptr),
      combigridMultiOperation(nullptr),
      combigridTensorOperation(combigridTensorOperation),
      functionBases(0),
      numGridPoints(0) {
  // create vector of function bases
  for (size_t idim = 0; idim < numDims; idim++) {
    functionBases.push_back(functionBasis);
  }

  initializeLinearTransformation();
  initializeTensorOperation(combigridTensorOperation->getPointHierarchies(),
                            combigridTensorOperation->getStorage(),
                            combigridTensorOperation->getLevelManager());
}

PolynomialStochasticCollocation::PolynomialStochasticCollocation(
    std::shared_ptr<sgpp::combigrid::CombigridOperation> combigridOperation,
    std::vector<std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D>>& functionBases)
    : numDims(combigridOperation->numDims()),
      combigridOperation(combigridOperation),
      combigridMultiOperation(nullptr),
      combigridTensorOperation(nullptr),
      functionBases(functionBases),
      numGridPoints(0) {
  // make sure that the number of dimensions match
  if (numDims != functionBases.size()) {
    throw sgpp::base::application_exception(
        "number of basis function do not match with the number of dimensions of the operation");
  }

  initializeLinearTransformation();
  initializeTensorOperation(combigridOperation->getPointHierarchies(),
                            combigridOperation->getStorage(),
                            combigridOperation->getLevelManager());
}

PolynomialStochasticCollocation::PolynomialStochasticCollocation(
    std::shared_ptr<sgpp::combigrid::CombigridMultiOperation> combigridMultiOperation,
    std::vector<std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D>>& functionBases)
    : numDims(combigridMultiOperation->numDims()),
      combigridOperation(nullptr),
      combigridMultiOperation(combigridMultiOperation),
      combigridTensorOperation(nullptr),
      functionBases(functionBases),
      numGridPoints(0) {
  // make sure that the number of dimensions match
  if (numDims != functionBases.size()) {
    throw sgpp::base::application_exception(
        "number of basis function do not match with the number of dimensions of the operation");
  }

  initializeLinearTransformation();
  initializeTensorOperation(combigridMultiOperation->getPointHierarchies(),
                            combigridMultiOperation->getStorage(),
                            combigridMultiOperation->getLevelManager());
}

PolynomialStochasticCollocation::PolynomialStochasticCollocation(
    std::shared_ptr<sgpp::combigrid::CombigridTensorOperation> combigridTensorOperation,
    std::vector<std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D>>& functionBases)
    : numDims(combigridTensorOperation->numDims()),
      combigridOperation(nullptr),
      combigridMultiOperation(nullptr),
      combigridTensorOperation(nullptr),
      functionBases(functionBases),
      numGridPoints(0) {
  // make sure that the number of dimensions match
  if (numDims != functionBases.size()) {
    throw sgpp::base::application_exception(
        "number of basis function do not match with the number of dimensions of the operation");
  }

  initializeLinearTransformation();
  initializeTensorOperation(combigridTensorOperation->getPointHierarchies(),
                            combigridTensorOperation->getStorage(),
                            combigridTensorOperation->getLevelManager());
}

#endif

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
}

#ifdef USE_DAKOTA
void PolynomialStochasticCollocation::initializeLinearTransformation() {
  sgpp::base::DataVector bounds(2 * numDims);
  for (size_t idim = 0; idim < numDims; idim++) {
    Pecos::RealRealPair bounds_idim = functionBases[idim]->getRandomVariable()->bounds();
    bounds[2 * idim] = bounds_idim.first;
    bounds[2 * idim + 1] = bounds_idim.second;
  }
  trans.initialize(bounds);
}
#endif

double PolynomialStochasticCollocation::quad(MultiIndex i, MultiIndex j) {
  double ans = 1.0;
  // Gauss quadrature in each dimension

  // performing Gauss-Legendre integration
  base::DataVector roots;
  base::DataVector quadratureweights;
  auto& quadRule = base::GaussLegendreQuadRule1D::getInstance();

  for (size_t idim = 0; idim < i.size(); idim++) {
    size_t degree_i = i[idim], degree_j = j[idim];
    auto functionBasis = functionBases[idim];
    size_t numGaussPoints = (degree_i + degree_j + 3) / 2;
    //    if (functionBasis->getConfiguration().polyParameters.type_ !=
    //        OrthogonalPolynomialBasisType::LEGENDRE) {
    //      numGaussPoints = std::min(quadRule.getMaxSupportedLevel(), numGaussPoints);
    //    }

    // do iterative 1d quadrature
    double sum_idim = 0.0;
    double sum_idim_old = 0.0;
    double err = 1e14;
    size_t iteration = 0;
    std::cout << degree_i << ", " << degree_j << std::endl;
    while (err > 1e-13 && numGaussPoints < quadRule.getMaxSupportedLevel()) {
      quadRule.getLevelPointsAndWeightsNormalized(numGaussPoints, roots, quadratureweights);

      sum_idim = 0.0;
      for (size_t i = 0; i < roots.getSize(); ++i) {
        double x_unit = roots[i], w = quadratureweights[i];
        double x_prob = trans.unitToProbabilistic(x_unit, idim);
        sum_idim += w * legendreBasis->evaluate(degree_i, x_unit) *
                    legendreBasis->evaluate(degree_j, x_unit) * functionBasis->pdf(x_prob);
      }

      sum_idim *= trans.vol(idim);
      if (iteration > 0) {
        err = std::abs(sum_idim_old - sum_idim);
        if (std::abs(sum_idim) > 1e-13) {
          err = err / std::abs(sum_idim);
        }
        std::cout << "  d:" << idim << "i:" << iteration << "; err=" << (err * std::abs(sum_idim))
                  << "; relerr=" << err << std::endl;
      }
      sum_idim_old = sum_idim;
      numGaussPoints += 1;
      iteration += 1;
    }

    ans *= sum_idim;
  }
  return ans;
}

double PolynomialStochasticCollocation::quad(MultiIndex i) {
  double ans = 1.0;
  // Gauss quadrature in each dimension

  // performing Gauss-Legendre integration
  base::DataVector roots;
  base::DataVector quadratureweights;
  auto& quadRule = base::GaussLegendreQuadRule1D::getInstance();

  for (size_t idim = 0; idim < i.size(); idim++) {
    size_t degree_i = i[idim];
    size_t numGaussPoints = (degree_i + 2) / 2;
    auto functionBasis = functionBases[idim];

    // do iterative 1d quadrature
    double sum_idim = 0.0;
    double sum_idim_old = 0.0;
    double err = 1e14;
    size_t iteration = 0;
    while (err > 1e-13 && numGaussPoints < quadRule.getMaxSupportedLevel()) {
      quadRule.getLevelPointsAndWeightsNormalized(numGaussPoints, roots, quadratureweights);

      sum_idim = 0.0;
      for (size_t i = 0; i < roots.getSize(); ++i) {
        double x_unit = roots[i], w = quadratureweights[i];
        double x_prob = trans.unitToProbabilistic(x_unit, idim);
        sum_idim += w * legendreBasis->evaluate(degree_i, x_unit) * functionBasis->pdf(x_prob);
      }

      sum_idim *= trans.vol(idim);
      if (iteration > 0) {
        err = std::abs(sum_idim_old - sum_idim);
        if (std::abs(sum_idim) > 1e-13) {
          err = err / std::abs(sum_idim);
        }
        std::cout << "  d:" << idim << "i:" << iteration << "; err=" << (err * std::abs(sum_idim))
                  << "; relerr=" << err << std::endl;
      }
      sum_idim_old = sum_idim;
      numGaussPoints += 1;
      iteration += 1;
    }

    ans *= sum_idim;
  }
  return ans;
}

double PolynomialStochasticCollocation::mean() {
  if (numGridPoints < combigridTensorOperation->numGridPoints()) {
    expansionCoefficients = combigridTensorOperation->getResult();
    numGridPoints = combigridTensorOperation->numGridPoints();
  }

  // compute mass matrix and corresponding coefficient vector
  auto multiIndices_i = expansionCoefficients.getValues();
  auto it_i = multiIndices_i->getStoredDataIterator();
  sgpp::base::DataVector m(numGridPoints);
  sgpp::base::DataVector coeffs(numGridPoints);

  size_t i = 0;
  while (it_i->isValid() && i < numGridPoints) {
    auto ix = it_i->getMultiIndex();
    m[i] = quad(ix);
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

double PolynomialStochasticCollocation::variance() {
  if (numGridPoints < combigridTensorOperation->numGridPoints()) {
    expansionCoefficients = combigridTensorOperation->getResult();
    numGridPoints = combigridTensorOperation->numGridPoints();
  }

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
        if (innerProducts.find(kx) != innerProducts.end()) {
          innerProduct = innerProducts[kx];
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

std::shared_ptr<sgpp::combigrid::CombigridTensorOperation>
PolynomialStochasticCollocation::getCombigridTensorOperation() {
  return combigridTensorOperation;
}

void PolynomialStochasticCollocation::joinMultiIndices(MultiIndex& ix, MultiIndex& jx,
                                                       MultiIndex& kx) {
  kx.insert(kx.end(), ix.begin(), ix.end());
  kx.insert(kx.end(), jx.begin(), jx.end());
}

} /* namespace combigrid */
} /* namespace sgpp */
