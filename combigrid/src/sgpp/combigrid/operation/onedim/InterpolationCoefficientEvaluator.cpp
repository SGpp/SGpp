// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/combigrid/functions/OrthogonalPolynomialBasis1D.hpp>
#include <sgpp/combigrid/operation/onedim/InterpolationCoefficientEvaluator.hpp>

#include <algorithm>
#include <iostream>
#include <vector>

#ifdef USE_EIGEN
#include <eigen3/Eigen/Dense>
#endif

namespace sgpp {
namespace combigrid {

InterpolationCoefficientEvaluator::InterpolationCoefficientEvaluator()
    : basisValues(1, FloatTensorVector(1)), basisCoefficients() {
  evalConfig.type = CombiEvaluatorTypes::Tensor_PolynomialInterpolation;
  sgpp::combigrid::OrthogonalPolynomialBasis1DConfiguration config;
  config.polyParameters.type_ = sgpp::combigrid::OrthogonalPolynomialBasisType::LEGENDRE;
  functionBasis = std::make_shared<OrthogonalPolynomialBasis1D>(config);
}

InterpolationCoefficientEvaluator::InterpolationCoefficientEvaluator(
    std::shared_ptr<AbstractInfiniteFunctionBasis1D> functionBasis)
    : basisValues(1, FloatTensorVector(1)), basisCoefficients(), functionBasis(functionBasis) {
  evalConfig.type = CombiEvaluatorTypes::Tensor_PolynomialInterpolation;
  evalConfig.functionBasis = functionBasis;
  if (functionBasis == nullptr) {
    throw sgpp::base::application_exception(
        "InterpolationCoefficientEvaluator: the basis function is not defined");
  }
}

InterpolationCoefficientEvaluator::InterpolationCoefficientEvaluator(
    const InterpolationCoefficientEvaluator& other)
    : basisValues(other.basisValues),
      basisCoefficients(other.basisCoefficients),
      functionBasis(other.functionBasis) {
  evalConfig.type = CombiEvaluatorTypes::Tensor_PolynomialInterpolation;
  evalConfig.functionBasis = other.functionBasis;
}

InterpolationCoefficientEvaluator::~InterpolationCoefficientEvaluator() {}

void InterpolationCoefficientEvaluator::setGridPoints(const std::vector<double>& xValues) {
#ifdef USE_EIGEN
  size_t n = xValues.size();
  Eigen::MatrixXd mat(n, n);

  // compute the interpolation matrix
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      mat(i, j) = functionBasis->evaluate(j, xValues[i]);
    }
  }

  // use Christoffel preconditioner: scale length of each row entry to 1 to reduce
  // the condition number of the interpolation matrix
  Eigen::VectorXd invMaxNorms = mat.rowwise().norm();
  for (size_t i = 0; i < n; ++i) {
    invMaxNorms(i) = 1. / invMaxNorms(i);
  }
  mat = mat.array().colwise() * invMaxNorms.array();

  // invert the preconditioned matrix
  Eigen::MatrixXd invertedScaledMat = mat.fullPivHouseholderQr().inverse();

  // undo the preconditioning to the inverse
  basisValues = std::vector<FloatTensorVector>(n, FloatTensorVector(1));
  for (size_t j = 0; j < n; ++j) {
    for (size_t i = 0; i < n; ++i) {
      basisValues[j].at(MultiIndex{i}) = invertedScaledMat(i, j) * invMaxNorms(j);
    }
  }

#else
  throw sgpp::base::application_exception("need Eigen to use the PCE transformation.");
#endif
}

void InterpolationCoefficientEvaluator::setBasisCoefficientsAtGridPoints(
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
