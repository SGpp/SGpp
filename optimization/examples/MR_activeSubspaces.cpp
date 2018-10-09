// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/optimization/function/scalar/WrapperScalarFunction.hpp>
#include "../src/sgpp/optimization/activeSubspaces/ASMatrixNakBsplineBoundary.hpp"
#include "../src/sgpp/optimization/activeSubspaces/ASResponseSurfaceNakBsplineModified.hpp"

double f(sgpp::base::DataVector v) { return exp(0.7 * v[0] + 0.3 * v[1]); }

int main() {
  size_t numDim = 2;
  size_t degree = 3;
  size_t level = 3;
  // active subspace specifier
  size_t n = 1;
  // set k to one more than the desired reduced dimensionality
  size_t k = n + 1;
  // choose an oversampling factor alpha between 2 and 10 (the larger, the more evaluations, the
  // better the approximation of the eigenvalues)
  size_t alpha = 3;
  // the number of samples ensuring that the first k eigenvalues are accurate enough
  size_t numMCPoints = alpha * k * static_cast<size_t>(ceil(log(static_cast<double>(numDim))));
  numMCPoints = 2000;

  sgpp::optimization::WrapperScalarFunction objectiveFunc(numDim, f);
  sgpp::optimization::ASMatrixNakBsplineBoundary ASM(objectiveFunc, degree);
  ASM.buildRegularInterpolant(level);
  ASM.createMatrix(numMCPoints);
  ASM.evDecomposition();

  Eigen::VectorXd eigenvalues = ASM.getEigenvalues();
  std::cout << "EVal:\n" << eigenvalues << std::endl;
  Eigen::MatrixXd eigenvectors = ASM.getEigenvectors();
  std::cout << "EVec:\n" << eigenvectors << std::endl;
  Eigen::MatrixXd W1 = ASM.getTransformationMatrix(n);
  std::cout << "W1:\n" << W1 << std::endl;

  //----------------------------
  Eigen::MatrixXd M = ASM.getMatrix();
  std::cout << "MC:\n" << M << std::endl;

  ASM.createMatrixGauss();
  M = ASM.getMatrix();
  std::cout << "Gauss:\n" << M << std::endl;
  //----------------------------

  sgpp::base::DataMatrix evaluationPoints = ASM.getEvaluationPoints();
  sgpp::base::DataVector functionValues = ASM.getFunctionValues();

  //  sgpp::optimization::ASResponseSurfaceNakBspline responseSurf(W1);
  //  size_t responseLevel = 5;
  //  responseSurf.createRegularSurfaceFromDetectionPoints(evaluationPoints, functionValues,
  //                                                       responseLevel);
  //
  //  sgpp::base::DataVector v(numDim, 0.5);
  //  std::cout << "f(v)  = " << f(v) << std::endl;
  //  std::cout << "rI(v) = " << responseSurf.eval(v) << std::endl;

  return 0;
}
