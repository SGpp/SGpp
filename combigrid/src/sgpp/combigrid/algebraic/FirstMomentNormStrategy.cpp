// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/algebraic/FirstMomentNormStrategy.hpp>
#include <sgpp/combigrid/integration/GaussLegendreQuadrature.hpp>
#include <sgpp/combigrid/functions/OrthogonalPolynomialBasis1D.hpp>
#include <sgpp/combigrid/functions/WeightFunctionsCollection.hpp>
#include <sgpp/combigrid/GeneralFunction.hpp>

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/combigrid/functions/OrthogonalBasisFunctionsCollection.hpp>

#include <vector>
#include <map>

namespace sgpp {
namespace combigrid {

FirstMomentNormStrategy::FirstMomentNormStrategy(
    std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D> basisFunction, size_t numDims,
    sgpp::combigrid::SingleFunction weightFunction, bool isOrthogonal,
    sgpp::base::DataVector const& bounds)
    : isOrthogonal(isOrthogonal),
      bounds(bounds),
      basisFunctions(numDims, basisFunction),
      weightFunctions(numDims, weightFunction) {
  initializeBounds();
}

FirstMomentNormStrategy::FirstMomentNormStrategy(
    std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D> basisFunction,
    sgpp::combigrid::WeightFunctionsCollection& weightFunctions, bool isOrthogonal,
    sgpp::base::DataVector const& bounds)
    : isOrthogonal(isOrthogonal),
      bounds(bounds),
      basisFunctions(weightFunctions.size(), basisFunction),
      weightFunctions(weightFunctions) {
  initializeBounds();
}

FirstMomentNormStrategy::FirstMomentNormStrategy(
    sgpp::combigrid::OrthogonalBasisFunctionsCollection& basisFunctions,
    sgpp::combigrid::WeightFunctionsCollection& weightFunctions, bool isOrthogonal,
    sgpp::base::DataVector const& bounds)
    : isOrthogonal(isOrthogonal),
      bounds(bounds),
      basisFunctions(basisFunctions),
      weightFunctions(weightFunctions) {
  initializeBounds();
}

FirstMomentNormStrategy::~FirstMomentNormStrategy() {}

double FirstMomentNormStrategy::quad(MultiIndex i, GaussLegendreQuadrature& quadRule) {
  double ans = 1.0;

  // Gauss quadrature in each dimension
  for (size_t idim = 0; idim < i.size(); idim++) {
    size_t degree_i = i[idim];
    auto basisFunction = basisFunctions[idim];
    auto weightFunction = weightFunctions[idim];
    size_t incrementQuadraturePoints = basisFunction->numAdditionalQuadraturePoints();
    size_t numGaussPoints = (degree_i + 2) / 2;

    auto func = [&basisFunction, &degree_i, &weightFunction](double x_unit, double x_prob) {
      return basisFunction->evaluate(degree_i, x_unit) * weightFunction(x_unit);
    };

    double a = bounds[2 * idim], b = bounds[2 * idim + 1];
    if (incrementQuadraturePoints == 0) {
      quadRule.initialize(numGaussPoints);
      ans *= GaussLegendreQuadrature(numGaussPoints).evaluate(func, a, b);
    } else {
      ans *= quadRule.evaluate_iteratively(func, a, b, numGaussPoints + incrementQuadraturePoints,
                                           incrementQuadraturePoints, 1e-14);
    }
  }
  return ans;
}

double FirstMomentNormStrategy::computeMean(FloatTensorVector& vector) {
  auto it_i = vector.getValues()->getStoredDataIterator();

  // initialize Gauss quadrature
  GaussLegendreQuadrature quadRule(100);

  double ans = 0.0;
  while (it_i->isValid()) {
    auto ix = it_i->getMultiIndex();

    double coeff = it_i->value().value();

    double quadValue = 0.0;
    auto it = lookupTable.find(ix);
    if (it != lookupTable.end()) {
      quadValue = it->second;
    } else {
      quadValue = quad(ix, quadRule);
      lookupTable[ix] = quadValue;
    }

    ans += quadValue * coeff;

    it_i->moveToNext();
  }

  return ans;
}

double FirstMomentNormStrategy::norm(FloatTensorVector& vector) {
  if (isOrthogonal) {
    auto values = vector.getValues();
    MultiIndex ix(values->getNumDimensions(), 0);
    return values->get(ix).value();
  } else {
    return computeMean(vector);
  }
}

void FirstMomentNormStrategy::initializeBounds() {
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
          "FirstMomentNormStrategy::initializeBounds - not enough arguments for bounds "
          "specified");
    }
  }
}

} /* namespace combigrid */
} /* namespace sgpp */
