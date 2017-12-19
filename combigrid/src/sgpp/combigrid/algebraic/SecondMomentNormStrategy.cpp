// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/algebraic/SecondMomentNormStrategy.hpp>
#include <sgpp/combigrid/integration/GaussLegendreQuadrature.hpp>
#include <sgpp/combigrid/functions/OrthogonalPolynomialBasis1D.hpp>

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <vector>

namespace sgpp {
namespace combigrid {

SecondMomentNormStrategy::SecondMomentNormStrategy(
    std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D> basisFunction, size_t numDims,
    bool isOrthogonal)
    : isOrthogonal(isOrthogonal), basisFunctions(numDims, basisFunction) {}
SecondMomentNormStrategy::SecondMomentNormStrategy(
    std::vector<std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D>>& basisFunctions,
    bool isOrthogonal)
    : isOrthogonal(isOrthogonal), basisFunctions(basisFunctions) {}
SecondMomentNormStrategy::~SecondMomentNormStrategy() {}

double SecondMomentNormStrategy::quad(MultiIndex i, MultiIndex j) {
  double ans = 1.0;

  // Gauss quadrature in each dimension
  for (size_t idim = 0; idim < i.size(); idim++) {
    size_t degree_i = i[idim];
    size_t degree_j = j[idim];
    auto functionBasis = basisFunctions[idim];
    size_t numGaussPoints = (degree_i + degree_j + 3) / 2;

    auto func = [functionBasis, &degree_i, &degree_j](double x_unit, double x_prob) {
      return functionBasis->evaluate(degree_i, x_unit) * functionBasis->evaluate(degree_j, x_unit) *
             functionBasis->pdf(x_prob);
    };

    ans *= GaussLegendreQuadrature::evaluate_iteratively(
        func, functionBasis->lowerBound(), functionBasis->upperBound(), numGaussPoints,
        functionBasis->numAdditionalQuadraturePoints());
  }
  return ans;
}

double SecondMomentNormStrategy::computeSecondMoment(FloatTensorVector& vector) {
  // compute mass matrix and corresponding coefficient vector
  auto multiIndices_i = vector.getValues();
  auto multiIndices_j = vector.getValues();
  sgpp::base::DataMatrix M(0, 0);
  sgpp::base::DataVector coeffs;

  auto it_i = multiIndices_i->getStoredDataIterator();
  size_t i = 0;
  while (it_i->isValid()) {
    MultiIndex ix = it_i->getMultiIndex();

    auto it_j = multiIndices_j->getStoredDataIterator();
    size_t j = 0;
    sgpp::base::DataVector row;
    while (it_j->isValid()) {
      double innerProduct = 0.0;
      // exploit symmetry
      if (j >= i) {
        MultiIndex jx = it_j->getMultiIndex();

        // compute the inner product and store it
        innerProduct = quad(ix, jx);
      } else {
        innerProduct = M.get(j, i);
      }
      row.push_back(innerProduct);

      it_j->moveToNext();
      j += 1;
    }
    M.appendRow(row);

    it_i->moveToNext();
    i += 1;
  }

  // compute second moment: coeffs^T M coeffs
  sgpp::base::DataVector result(coeffs.getSize());
  M.mult(coeffs, result);

  double ans = 0.0;
  for (size_t i = 0; i < M.getNrows(); i++) {
    ans += coeffs[i] * result[i];
  }

  return ans;
}

} /* namespace combigrid */
} /* namespace sgpp */
