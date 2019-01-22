// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

// #ifdef USE_EIGEN

#include <sgpp/optimization/activeSubspaces/ASMatrixBsplineData.hpp>

namespace sgpp {
namespace optimization {

void ASMatrixBsplineData::buildRegularInterpolant(size_t level) {
  grid->getGenerator().regular(level);
  this->calculateInterpolationCoefficients();
}

void ASMatrixBsplineData::buildAdaptiveInterpolant(size_t maxNumGridPoints, size_t initialLevel,
                                                   size_t refinementsNum) {
  grid->getGenerator().regular(initialLevel);
  this->calculateInterpolationCoefficients();
  while (grid->getSize() < maxNumGridPoints) {
    this->refineSurplusAdaptive(refinementsNum);
  }
}

void ASMatrixBsplineData::calculateInterpolationCoefficients() {
  // regressionMatrix(i,j) = b_j (x_i)
  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  Eigen::MatrixXd regressionMatrix(evaluationPoints.getNcols(), gridStorage.getSize());
  for (size_t i = 0; i < evaluationPoints.getNcols(); i++) {
    sgpp::base::DataVector point(evaluationPoints.getNrows());
    evaluationPoints.getColumn(i, point);
    Eigen::VectorXd y_i = DataVectorToEigen(point);
    for (size_t j = 0; j < gridStorage.getSize(); j++) {
      double basisEval = 1;
      for (size_t t = 0; t < gridStorage.getDimension(); t++) {
        double basisEval1D =
            basis->eval(gridStorage.getPointLevel(j, t), gridStorage.getPointIndex(j, t), y_i(t));
        if (basisEval1D == 0) {
          basisEval = 0;
          break;
        } else {
          basisEval *= basisEval1D;
        }
      }
      regressionMatrix(i, j) = basisEval;
    }
  }
  Eigen::VectorXd functionValues_Eigen = DataVectorToEigen(functionValues);

  // three different ways to solve least squares with Eigen
  // https://eigen.tuxfamily.org/dox/group__LeastSquares.html

  Eigen::setNbThreads(4);
  /* 1. SVD, slowest but most accurate*/
  //  Eigen::VectorXd alpha_Eigen = regressionMatrix.bdcSvd(Eigen::ComputeThinU |
  //  Eigen::ComputeThinV)
  //                                    .solve(functionValues_Eigen);
  /* 2. QR medium speed, medium accuracy*/
  //  Eigen::VectorXd alpha_Eigen =
  //  regressionMatrix.colPivHouseholderQr().solve(functionValues_Eigen);
  /* 3. normal equations, fastest but least accurate. Terrible for ill conditioned matrices*/
  Eigen::VectorXd alpha_Eigen = (regressionMatrix.transpose() * regressionMatrix)
                                    .ldlt()
                                    .solve(regressionMatrix.transpose() * functionValues_Eigen);

  coefficients = EigenToDataVector(alpha_Eigen);
}

double ASMatrixBsplineData::l2InterpolationError(Eigen::MatrixXd errorPoints,
                                                 Eigen::VectorXd values) {
  sgpp::optimization::InterpolantScalarFunction interpolant(*grid, coefficients);
  Eigen::VectorXd evaluations(errorPoints.cols());
  for (unsigned int i = 0; i < errorPoints.cols(); i++) {
    evaluations(i) = interpolant.eval(EigenToDataVector(errorPoints.col(i)));
  }
  return (evaluations - values).norm();
}

}  // namespace optimization
}  // namespace sgpp

// #endif /* USE_EIGEN */
