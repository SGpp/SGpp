// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/algebraic/FirstMomentNormStrategy.hpp>
#include <sgpp/combigrid/integration/GaussLegendreQuadrature.hpp>
#include <sgpp/combigrid/functions/OrthogonalPolynomialBasis1D.hpp>

#include <vector>

namespace sgpp {
namespace combigrid {

FirstMomentNormStrategy::FirstMomentNormStrategy(
    std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D> basisFunction, size_t numDims,
    bool isOrthogonal)
    : isOrthogonal(isOrthogonal), basisFunctions(numDims, basisFunction) {}
FirstMomentNormStrategy::FirstMomentNormStrategy(
    std::vector<std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D>>& basisFunctions,
    bool isOrthogonal)
    : isOrthogonal(isOrthogonal), basisFunctions(basisFunctions) {}

FirstMomentNormStrategy::~FirstMomentNormStrategy() {}

double FirstMomentNormStrategy::quad(MultiIndex i) {
  double ans = 1.0;

  // Gauss quadrature in each dimension
  for (size_t idim = 0; idim < i.size(); idim++) {
    size_t degree_i = i[idim];
    auto functionBasis = basisFunctions[idim];
    size_t incrementQuadraturePoints = functionBasis->numAdditionalQuadraturePoints();
    size_t numGaussPoints = (degree_i + 2) / 2;

    auto func = [&functionBasis, &degree_i, &idim, this](double x_unit, double x_prob) {
      return functionBasis->evaluate(degree_i, x_unit) * functionBasis->pdf(x_prob);
    };

    ans *= GaussLegendreQuadrature::evaluate_iteratively(func, functionBasis->lowerBound(),
                                                         functionBasis->upperBound(),
                                                         numGaussPoints, incrementQuadraturePoints);
  }
  return ans;
}

double FirstMomentNormStrategy::computeMean(FloatTensorVector& vector) {
  // compute mass matrix and corresponding coefficient vector
  auto multiIndices_i = vector.getValues();
  auto it_i = multiIndices_i->getStoredDataIterator();
  sgpp::base::DataVector m;
  sgpp::base::DataVector coeffs;

  while (it_i->isValid()) {
    auto ix = it_i->getMultiIndex();

    double value = quad(ix);
    m.push_back(value);
    coeffs.push_back(it_i->value().value());
    it_i->moveToNext();
  }

  // compute mean: m^T coeffs
  double ans = 0.0;
  for (size_t i = 0; i < m.getSize(); i++) {
    ans += m[i] * coeffs[i];
  }

  return ans;
}

} /* namespace combigrid */
} /* namespace sgpp */
