// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

// #ifdef USE_EIGEN

#include <sgpp/datadriven/activeSubspaces/ASMatrixBsplineData.hpp>

namespace sgpp {
namespace datadriven {

void ASMatrixBsplineData::buildRegularInterpolant(size_t level) {
  grid->getGenerator().regular(level);
  this->calculateCoefficients();
}

void ASMatrixBsplineData::buildAdaptiveInterpolant(size_t maxNumGridPoints, size_t initialLevel,
                                                   size_t refinementsNum) {
  grid->getGenerator().regular(initialLevel);
  this->calculateCoefficients();
  while (grid->getSize() < maxNumGridPoints) {
    this->refineSurplusAdaptive(refinementsNum);
  }
}

void ASMatrixBsplineData::calculateCoefficients() {
  double mse = 0;
  sgpp::base::DataVector errorPerBasis;
  coefficients = EigenRegression(grid, degree, DataMatrixToEigen(evaluationPoints), functionValues,
                                 mse, errorPerBasis);
}

double ASMatrixBsplineData::l2InterpolationError(Eigen::MatrixXd errorPoints,
                                                 Eigen::VectorXd values) {
  sgpp::base::InterpolantScalarFunction interpolant(*grid, coefficients);
  Eigen::VectorXd evaluations(errorPoints.cols());
  for (unsigned int i = 0; i < errorPoints.cols(); i++) {
    evaluations(i) = interpolant.eval(EigenToDataVector(errorPoints.col(i)));
  }
  return (evaluations - values).norm();
}

}  // namespace datadriven
}  // namespace sgpp

// #endif /* USE_EIGEN */
