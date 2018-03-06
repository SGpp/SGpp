// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/algebraic/SecondMomentNormStrategy.hpp>
#include <sgpp/combigrid/integration/GaussLegendreQuadrature.hpp>
#include <sgpp/combigrid/functions/OrthogonalPolynomialBasis1D.hpp>
#include <sgpp/combigrid/functions/WeightFunctionsCollection.hpp>
#include <sgpp/combigrid/GeneralFunction.hpp>

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/combigrid/functions/OrthogonalBasisFunctionsCollection.hpp>

#include <vector>

namespace sgpp {
namespace combigrid {

SecondMomentNormStrategy::SecondMomentNormStrategy(
    std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D> basisFunction, size_t numDims,
    sgpp::combigrid::SingleFunction weightFunction, bool isOrthogonal,
    sgpp::base::DataVector const& bounds)
    : isOrthogonal(isOrthogonal),
      bounds(bounds),
      basisFunctions(numDims, basisFunction),
      weightFunctions(numDims, weightFunction) {
  initializeBounds();
}

SecondMomentNormStrategy::SecondMomentNormStrategy(
    std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D> basisFunction,
    sgpp::combigrid::WeightFunctionsCollection& weightFunctions, bool isOrthogonal,
    sgpp::base::DataVector const& bounds)
    : isOrthogonal(isOrthogonal),
      bounds(bounds),
      basisFunctions(weightFunctions.size(), basisFunction),
      weightFunctions(weightFunctions) {
  initializeBounds();
}

SecondMomentNormStrategy::SecondMomentNormStrategy(
    sgpp::combigrid::OrthogonalBasisFunctionsCollection& basisFunctions,
    sgpp::combigrid::WeightFunctionsCollection& weightFunctions, bool isOrthogonal,
    sgpp::base::DataVector const& bounds)
    : isOrthogonal(isOrthogonal),
      bounds(bounds),
      basisFunctions(basisFunctions),
      weightFunctions(weightFunctions) {
  initializeBounds();
}

SecondMomentNormStrategy::~SecondMomentNormStrategy() {}

double SecondMomentNormStrategy::quad(MultiIndex i, MultiIndex j,
                                      GaussLegendreQuadrature& quadRule) {
  double ans = 1.0;

  // Gauss quadrature in each dimension
  for (size_t idim = 0; idim < i.size(); idim++) {
    size_t degree_i = i[idim];
    size_t degree_j = j[idim];
    auto basisFunction = basisFunctions[idim];
    auto weightFunction = weightFunctions[idim];
    size_t incrementQuadraturePoints = basisFunction->numAdditionalQuadraturePoints();
    size_t numGaussPoints = (degree_i + degree_j + 3) / 2;

    auto func = [basisFunction, &degree_i, &degree_j, &weightFunction](double x_unit,
                                                                       double x_prob) {
      return basisFunction->evaluate(degree_i, x_unit) * basisFunction->evaluate(degree_j, x_unit) *
             weightFunction(x_unit);
    };

    double a = bounds[2 * idim], b = bounds[2 * idim + 1];
    if (incrementQuadraturePoints == 0) {
      quadRule.initialize(numGaussPoints);
      ans *= quadRule.evaluate(func, a, b);
    } else {
      ans *= quadRule.evaluate_iteratively(func, a, b, numGaussPoints + incrementQuadraturePoints,
                                           incrementQuadraturePoints, 1e-13);
    }
  }
  return ans;
}

double SecondMomentNormStrategy::computeSecondMoment(FloatTensorVector& vector) {
  // update the missing entries in the lookup table
  // compute mass matrix and corresponding coefficient vector
  auto multiIndices_i = vector.getValues();
  auto multiIndices_j = vector.getValues();

  // initialize Gauss quadrature
  GaussLegendreQuadrature quadRule(100);

  // count the number of tensor terms and compute coefficient vector
  auto it_counter = multiIndices_i->getStoredDataIterator();

  sgpp::base::DataVector coeffs;
  size_t numGridPoints = 0;
  while (it_counter->isValid()) {
    coeffs.push_back(it_counter->value().value());
    numGridPoints++;
    it_counter->moveToNext();
  }

  auto it_i = multiIndices_i->getStoredDataIterator();
  size_t i = 0;
  size_t numDims = multiIndices_i->getNumDimensions();
  MultiIndex kx(2 * numDims);

  sgpp::base::DataMatrix M(numGridPoints, numGridPoints);
  sgpp::base::DataVector result(numGridPoints);
  while (it_i->isValid()) {
    MultiIndex ix = it_i->getMultiIndex();

    auto it_j = multiIndices_j->getStoredDataIterator();
    size_t j = 0;
    while (it_j->isValid()) {
      // exploit symmetry
      if (j >= i) {
        MultiIndex jx = it_j->getMultiIndex();

        joinMultiIndices(ix, jx, kx);
        double innerProduct = 0.0;
        auto it = innerProducts.find(kx);
        if (it != innerProducts.end()) {
          innerProduct = it->second;
        } else {
          innerProduct = quad(ix, jx, quadRule);
          innerProducts[kx] = innerProduct;
        }

        // compute the matrix vector product
        result[i] += innerProduct * coeffs[j];
        if (j > i) {
          result[j] += innerProduct * coeffs[i];
        }
      }

      it_j->moveToNext();
      j += 1;
    }

    it_i->moveToNext();
    i += 1;
  }

  // compute the vector vector product
  double ans = 0.0;
  for (size_t i = 0; i < numGridPoints; i++) {
    ans += coeffs[i] * result[i];
  }

  return ans;
}

double SecondMomentNormStrategy::norm(FloatTensorVector& vector) {
  if (isOrthogonal) {
    double sum = 0.0;
    for (auto it = vector.getValues()->getStoredDataIterator(); it->isValid(); it->moveToNext()) {
      double coeff = it->value().value();
      sum += coeff * coeff;
    }
    return sum;
  } else {
    return computeSecondMoment(vector);
  }
}

void SecondMomentNormStrategy::initializeBounds() {
  size_t numDims = basisFunctions.size();
  if (bounds.size() == 0) {
    bounds.resize(2 * numDims);
    for (size_t idim = 0; idim < numDims; idim++) {
      bounds[2 * idim] = basisFunctions[idim]->lowerBound();
      bounds[2 * idim + 1] = basisFunctions[idim]->upperBound();
    }
  } else {
  }
}

void SecondMomentNormStrategy::joinMultiIndices(MultiIndex& ix, MultiIndex& jx, MultiIndex& kx) {
  for (size_t i = 0; i < ix.size(); i++) {
    if (ix[i] < jx[i]) {
      kx[2 * i] = ix[i];
      kx[2 * i + 1] = jx[i];
    } else {
      kx[2 * i] = jx[i];
      kx[2 * i + 1] = ix[i];
    }
  }
}

} /* namespace combigrid */
} /* namespace sgpp */
