// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/exception/generation_exception.hpp>
#include <sgpp/combigrid/functions/OrthogonalPolynomialBasis1D.hpp>
#include <sgpp/combigrid/operation/onedim/InterpolationCoefficientEvaluator.hpp>

#include <vector>

#ifdef USE_EIGEN
#include <eigen3/Eigen/Dense>
#endif

namespace sgpp {
namespace combigrid {

InterpolationCoefficientEvaluator::InterpolationCoefficientEvaluator()
    : basisValues(1, FloatTensorVector(1)), basisCoefficients() {
  sgpp::combigrid::OrthogonalPolynomialBasis1DConfiguration config;
  config.polyParameters.type_ = sgpp::combigrid::OrthogonalPolynomialBasisType::LEGENDRE;
  functionBasis = std::make_shared<OrthogonalPolynomialBasis1D>(config);
}

InterpolationCoefficientEvaluator::InterpolationCoefficientEvaluator(
    std::shared_ptr<AbstractInfiniteFunctionBasis1D> functionBasis)
    : basisValues(1, FloatTensorVector(1)), basisCoefficients(), functionBasis(functionBasis) {}

InterpolationCoefficientEvaluator::~InterpolationCoefficientEvaluator() {}

InterpolationCoefficientEvaluator::InterpolationCoefficientEvaluator(
    const InterpolationCoefficientEvaluator& other)
    : basisValues(other.basisValues),
      basisCoefficients(other.basisCoefficients),
      functionBasis(other.functionBasis) {}

void InterpolationCoefficientEvaluator::setGridPoints(const std::vector<double>& xValues) {
#ifdef USE_EIGEN
  size_t n = xValues.size();
  Eigen::MatrixXd mat(n, n);

  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      mat(i, j) = functionBasis->evaluate(j, xValues[i]);
    }
  }
  Eigen::MatrixXd invertedMatrix = mat.inverse();
  basisValues = std::vector<FloatTensorVector>(n, FloatTensorVector(1));
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      basisValues[j].at(MultiIndex{i}) = invertedMatrix(i, j);
    }
  }
#else
  throw new sgpp::base::generation_exception("need Eigen to use the PCE transformation.");
#endif
}

void InterpolationCoefficientEvaluator::setFunctionValuesAtGridPoints(
    std::vector<double>& functionValues) {
  basisCoefficients = functionValues;
}

std::shared_ptr<AbstractLinearEvaluator<FloatTensorVector> >
InterpolationCoefficientEvaluator::cloneLinear() {
  return std::make_shared<InterpolationCoefficientEvaluator>(*this);
}

bool InterpolationCoefficientEvaluator::needsOrderedPoints() { return false; }

bool InterpolationCoefficientEvaluator::needsParameter() { return false; }

void InterpolationCoefficientEvaluator::setParameter(const FloatTensorVector& param) {
  // do nothing because no parameter is needed
}

} /* namespace combigrid */
} /* namespace sgpp */
