// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/operation/onedim/BsplineInterpolationCoefficientEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/BSplineInterpolationEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/BSplineScalarProductEvaluator.hpp>

#include <sgpp/base/exception/generation_exception.hpp>

#include <vector>
#include <algorithm>
#include <iostream>

#ifdef USE_EIGEN
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Cholesky>
#endif

namespace sgpp {
namespace combigrid {

BsplineInterpolationCoefficientEvaluator::BsplineInterpolationCoefficientEvaluator()
    : basisValues(1, FloatTensorVector(1)), basisCoefficients() {}

BsplineInterpolationCoefficientEvaluator::BsplineInterpolationCoefficientEvaluator(
    const BsplineInterpolationCoefficientEvaluator& other)
    : basisValues(other.basisValues), basisCoefficients(other.basisCoefficients) {}

BsplineInterpolationCoefficientEvaluator::~BsplineInterpolationCoefficientEvaluator() {}

void BsplineInterpolationCoefficientEvaluator::setGridPoints(const std::vector<double>& xValues) {
#ifdef USE_EIGEN
  size_t n = xValues.size();

  // compute the Gramian matrix
  BSplineScalarProductEvaluator scalarProductEval;
  scalarProductEval.setGridPoints(xValues);
  std::vector<FloatArrayVector> values = scalarProductEval.getBasisValues();

  Eigen::MatrixXd G(n, n);
  for (size_t i = 0; i < n; ++i) {
    FloatArrayVector result = values[i];
    for (size_t j = 0; j < n; ++j) {
      G(i, j) = result.get(j).getValue();
    }
  }

  // compute the Cholesky decomposition
  Eigen::MatrixXd L(G.llt().matrixL());

  // invert the Cholesky matrix
  Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(n, n);
  Eigen::MatrixXd C = L.triangularView<Eigen::Lower>().solve(identity);

  //  std::cout << "----------------------------------" << std::endl;
  //  std::cout << "C" << std::endl;
  //  for (size_t i = 0; i < n; i++) {
  //    std::cout << "[";
  //    for (size_t j = 0; j < n; j++) {
  //      std::cout << C(i, j) << ", ";
  //    }
  //    std::cout << "]" << std::endl << " ";
  //  }
  //  std::cout << "]" << std::endl;
  //  std::cout << "----------------------------------" << std::endl;

  Eigen::MatrixXd mat(n, n);

  BSplineInterpolationEvaluator evaluator;
  evaluator.setGridPoints(xValues);

  // compute the interpolation matrix
  for (size_t i = 0; i < n; ++i) {
    evaluator.setParameter(FloatScalarVector(xValues[i]));
    auto result = evaluator.getBasisValues();
    for (size_t j = 0; j < n; ++j) {
      mat(i, j) = result[j].getValue();
    }
  }

  // compute the interpolation matrix for the orthogonal bspline basis
  mat = mat * C.transpose();

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
  throw new sgpp::base::generation_exception("need Eigen to use the PCE transformation.");
#endif
}

void BsplineInterpolationCoefficientEvaluator::setBasisCoefficientsAtGridPoints(
    std::vector<double>& functionValues) {
  basisCoefficients = functionValues;
}

std::shared_ptr<AbstractLinearEvaluator<FloatTensorVector> >
BsplineInterpolationCoefficientEvaluator::cloneLinear() {
  return std::make_shared<BsplineInterpolationCoefficientEvaluator>(*this);
}

bool BsplineInterpolationCoefficientEvaluator::needsOrderedPoints() { return false; }

bool BsplineInterpolationCoefficientEvaluator::needsParameter() { return false; }

void BsplineInterpolationCoefficientEvaluator::setParameter(const FloatTensorVector& param) {
  // do nothing because no parameter is needed
}

} /* namespace combigrid */
} /* namespace sgpp */
