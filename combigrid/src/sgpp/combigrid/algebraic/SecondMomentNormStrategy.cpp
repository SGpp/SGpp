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

double SecondMomentNormStrategy::quad(MultiIndex i, MultiIndex j) {
  double ans = 1.0;

  // Gauss quadrature in each dimension
  for (size_t idim = 0; idim < i.size(); idim++) {
    size_t degree_i = i[idim];
    size_t degree_j = j[idim];
    auto basisFunction = basisFunctions[idim];
    auto weightFunction = weightFunctions[idim];
    size_t numGaussPoints = (degree_i + degree_j + 3) / 2;

    auto func = [basisFunction, &degree_i, &degree_j, &weightFunction](double x_unit,
                                                                       double x_prob) {
      return basisFunction->evaluate(degree_i, x_unit) * basisFunction->evaluate(degree_j, x_unit) *
             weightFunction(x_prob);
    };

    double a = bounds[2 * idim], b = bounds[2 * idim + 1];
    ans *= GaussLegendreQuadrature::evaluate_iteratively(
        func, a, b, numGaussPoints, basisFunction->numAdditionalQuadraturePoints(), 1e-14);
  }
  return ans;
}

double SecondMomentNormStrategy::computeSecondMoment(FloatTensorVector& vector) {
  // compute mass matrix and corresponding coefficient vector
  auto multiIndices_i = vector.getValues();
  auto multiIndices_j = vector.getValues();

  // count the number of tensor terms and compute coefficient vector
  auto it_counter = multiIndices_i->getStoredDataIterator();
  size_t numGridPoints = 0;
  sgpp::base::DataVector coeffs;
  while (it_counter->isValid()) {
    coeffs.push_back(it_counter->value().value());
    numGridPoints++;
    it_counter->moveToNext();
  }

  sgpp::base::DataMatrix M(numGridPoints, numGridPoints);
  auto it_i = multiIndices_i->getStoredDataIterator();
  size_t i = 0;
  while (it_i->isValid()) {
    MultiIndex ix = it_i->getMultiIndex();

    auto it_j = multiIndices_j->getStoredDataIterator();
    size_t j = 0;
    while (it_j->isValid()) {
      double innerProduct = 0.0;
      // exploit symmetry
      if (j >= i) {
        MultiIndex jx = it_j->getMultiIndex();

        // compute the inner product and store it
        innerProduct = quad(ix, jx);

        // store the result
        M.set(i, j, innerProduct);
        M.set(j, i, innerProduct);
      }

      it_j->moveToNext();
      j += 1;
    }

    it_i->moveToNext();
    i += 1;
  }

  // compute second moment: coeffs^T M coeffs
  sgpp::base::DataVector result(numGridPoints);
  M.mult(coeffs, result);

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
    if (bounds.size() != 2 * numDims) {
      throw sgpp::base::application_exception(
          "SecondMomentNormStrategy::initializeBounds - not enough arguments for bounds "
          "specified");
    }
  }
}

} /* namespace combigrid */
} /* namespace sgpp */
