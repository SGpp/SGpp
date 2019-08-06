// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/operation/onedim/BSplineInterpolationCoefficientEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/BSplineInterpolationEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/BSplineScalarProductEvaluator.hpp>
#include <sgpp/combigrid/utils/BSplineRoutines.hpp>
#include <sgpp/combigrid/utils/CombigridBSplineBasis.hpp>

#include <sgpp/base/exception/application_exception.hpp>

#include <algorithm>
#include <functional>
#include <iostream>
#include <vector>

#ifdef USE_EIGEN
#include <eigen3/Eigen/Cholesky>
#include <eigen3/Eigen/Dense>
#endif

namespace sgpp {
namespace combigrid {

size_t getGridLevelForExpUniformBoundaryGrid(size_t numGridPoints) {
  if (numGridPoints > 1) {
    return sgpp::combigrid::log2(numGridPoints - 2);
  } else {
    return 0;
  }
}

BSplineInterpolationCoefficientEvaluator::BSplineInterpolationCoefficientEvaluator()
    : basisValues(1, FloatTensorVector(1)), basisCoefficients(), degree(3) {}

BSplineInterpolationCoefficientEvaluator::BSplineInterpolationCoefficientEvaluator(size_t degree)
    : basisValues(1, FloatTensorVector(1)), basisCoefficients(), degree(degree) {}

BSplineInterpolationCoefficientEvaluator::BSplineInterpolationCoefficientEvaluator(
    const BSplineInterpolationCoefficientEvaluator& other)
    : basisValues(other.basisValues),
      basisCoefficients(other.basisCoefficients),
      degree(other.degree) {}

BSplineInterpolationCoefficientEvaluator::~BSplineInterpolationCoefficientEvaluator() {}

void BSplineInterpolationCoefficientEvaluator::setGridPoints(const std::vector<double>& xValues) {
#ifdef USE_EIGEN
  size_t n = xValues.size();

  //  std::cout << "----------------------------------" << std::endl;
  //  std::cout << "x" << std::endl;
  //  std::cout << "[";
  //  for (size_t j = 0; j < n - 1; j++) {
  //    std::cout << xValues[j] << ", ";
  //  }
  //  std::cout << xValues[xValues.size() - 1] << "]" << std::endl;
  //  std::cout << "----------------------------------" << std::endl;

  // compute the Gramian matrix
  BSplineScalarProductEvaluator scalarProductEval(degree);
  scalarProductEval.setGridPoints(xValues);
  std::vector<FloatArrayVector> values = scalarProductEval.getBasisValues();

  Eigen::MatrixXd G(n, n);
  for (size_t i = 0; i < n; ++i) {
    FloatArrayVector result = values[i];
    for (size_t j = 0; j < n; ++j) {
      G(i, j) = result.get(j).getValue();
    }
  }

  //  std::cout << "----------------------------------" << std::endl;
  //  std::cout << "G" << std::endl;
  //  std::cout << "[";
  //  for (size_t i = 0; i < n; i++) {
  //    std::cout << "[";
  //    for (size_t j = 0; j < n; j++) {
  //      std::cout << std::setw(15) << std::setprecision(10) << G(i, j) << ", ";
  //    }
  //    std::cout << "]" << std::endl << " ";
  //  }
  //  std::cout << "]" << std::endl;
  //  std::cout << "----------------------------------" << std::endl;

  // compute the Cholesky decomposition
  Eigen::MatrixXd L(G.llt().matrixL());

  // invert the Cholesky matrix
  Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(n, n);
  Eigen::MatrixXd C = L.triangularView<Eigen::Lower>().solve(identity);

  //  std::cout << "----------------------------------" << std::endl;
  //  std::cout << "C" << std::endl;
  //  std::cout << "[";
  //  for (size_t i = 0; i < n; i++) {
  //    std::cout << "[";
  //    for (size_t j = 0; j < n; j++) {
  //      std::cout << std::setw(15) << std::setprecision(10) << C(i, j) << ", ";
  //    }
  //    std::cout << "]" << std::endl << " ";
  //  }
  //  std::cout << "]" << std::endl;
  //  std::cout << "----------------------------------" << std::endl;

  Eigen::MatrixXd mat(n, n);

  BSplineInterpolationEvaluator evaluator(degree);
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

  std::cout << "----------------------------------" << std::endl;
  std::cout << "M" << std::endl;
  std::cout << "[";
  for (size_t i = 0; i < n; i++) {
    std::cout << "[";
    for (size_t j = 0; j < n; j++) {
      std::cout << std::setw(15) << std::setprecision(10) << mat(i, j) << ", ";
    }
    std::cout << "]" << std::endl << " ";
  }
  std::cout << "]" << std::endl;
  std::cout << "----------------------------------" << std::endl;

  //  // use Christoffel preconditioner: scale length of each row entry to 1 to reduce
  //  // the condition number of the interpolation matrix
  //  Eigen::VectorXd invMaxNorms = mat.rowwise().norm();
  //  for (size_t i = 0; i < n; ++i) {
  //    invMaxNorms(i) = 1. / invMaxNorms(i);
  //  }
  //  mat = mat.array().colwise() * invMaxNorms.array();
  //
  //  // invert the preconditioned matrix
  //  Eigen::MatrixXd invertedScaledMat = mat.fullPivHouseholderQr().inverse();

  // undo the preconditioning to the inverse
  size_t level = getGridLevelForExpUniformBoundaryGrid(n);
  size_t offset = getUniqueIndex(level, 0);
  basisValues = std::vector<FloatTensorVector>(n, FloatTensorVector(1));
  for (size_t j = 0; j < n; ++j) {
    for (size_t i = 0; i < n; ++i) {
      std::cout << "(" << (offset + i) << ") -> " << mat(i, j) << std::endl;
      basisValues[j].at(MultiIndex{offset + i}) = mat(i, j);
    }
  }

#else
  throw sgpp::base::application_exception("need Eigen to use the PCE transformation.");
#endif
}

void BSplineInterpolationCoefficientEvaluator::setBasisCoefficientsAtGridPoints(
    std::vector<double>& functionValues) {
  basisCoefficients = functionValues;
}

std::shared_ptr<AbstractLinearEvaluator<FloatTensorVector> >
BSplineInterpolationCoefficientEvaluator::cloneLinear() {
  return std::make_shared<BSplineInterpolationCoefficientEvaluator>(*this);
}

bool BSplineInterpolationCoefficientEvaluator::needsOrderedPoints() { return true; }

bool BSplineInterpolationCoefficientEvaluator::needsParameter() { return false; }

void BSplineInterpolationCoefficientEvaluator::setParameter(const FloatTensorVector& param) {
  // do nothing because no parameter is needed
}

} /* namespace combigrid */
} /* namespace sgpp */
