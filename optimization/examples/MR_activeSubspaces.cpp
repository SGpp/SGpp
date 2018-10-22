// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/optimization/activeSubspaces/ASMatrixGradientMC.hpp>
#include <sgpp/optimization/activeSubspaces/ASMatrixNakBspline.hpp>
#include <sgpp/optimization/activeSubspaces/ASResponseSurfaceNakBspline.hpp>
#include <sgpp/optimization/function/scalar/WrapperScalarFunction.hpp>
#include <sgpp/optimization/tools/Printer.hpp>

double f(sgpp::base::DataVector v) {
  //  return v[0] + 2 * v[1];
  return v[0] * v[0] + v[0] * v[1];
  //  return exp(0.7 * v[0] + 0.3 * v[1]);
}
double df(sgpp::base::DataVector v, sgpp::base::DataVector& gradient) {
  gradient.resizeZero(2);
  //  gradient[0] = 1;
  //  gradient[1] = 2;
  //  return v[0] + 2 * v[1];
  gradient[0] = 2 * v[0] + v[1];
  gradient[1] = v[0];
  return v[0] * v[0] + v[0] * v[1];
  //  gradient[0] = 0.7 * exp(0.7 * v[0] + 0.3 * v[1]);
  //  gradient[1] = 0.3 * exp(0.7 * v[0] + 0.3 * v[1]);
  //  return exp(0.7 * v[0] + 0.3 * v[1]);
}

int main() {
  sgpp::optimization::Printer::getInstance().setVerbosity(-1);
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
  sgpp::base::GridType gridType = sgpp::base::GridType::NakBsplineBoundary;
  sgpp::optimization::ASMatrixNakBspline ASM(objectiveFunc, gridType, degree);
  ASM.buildRegularInterpolant(level);
  ASM.createMatrixMonteCarlo(numMCPoints);
  ASM.evDecompositionForSymmetricMatrices();

  Eigen::VectorXd eigenvalues = ASM.getEigenvalues();
  std::cout << "EVal:\n" << eigenvalues << std::endl;
  Eigen::MatrixXd eigenvectors = ASM.getEigenvectors();
  //  std::cout << "EVec:\n" << eigenvectors << std::endl;
  Eigen::MatrixXd W1 = ASM.getTransformationMatrix(n);
  //  std::cout << "W1:\n" << W1 << std::endl;

  //----------------------------
  Eigen::MatrixXd M = ASM.getMatrix();
  std::cout << "MC (B-spline):\n" << M << std::endl;

  ASM.createMatrixGauss();
  M = ASM.getMatrix();
  std::cout << "Gauss (B-spline):\n" << M << std::endl;

  sgpp::optimization::WrapperScalarFunctionGradient objectiveFuncGradient(numDim, df);
  sgpp::optimization::ASMatrixGradientMC ASMGradient(objectiveFunc, objectiveFuncGradient);
  ASMGradient.createMatrixMonteCarlo(numMCPoints);
  M = ASMGradient.getMatrix();
  std::cout << "MC (Gradient):\n" << M << std::endl;
  //----------------------------

  sgpp::base::DataMatrix evaluationPoints = ASM.getEvaluationPoints();
  sgpp::base::DataVector functionValues = ASM.getFunctionValues();

  sgpp::optimization::ASResponseSurfaceNakBspline responseSurf(W1, gridType, degree);
  size_t responseLevel = 5;
  //  responseSurf.createRegularSurfaceFromDetectionPoints(evaluationPoints, functionValues,
  //                                                       responseLevel);
  responseSurf.createRegularResponseSurfaceWithPseudoInverse(responseLevel, objectiveFunc);

  sgpp::base::DataVector v(numDim, 0.3371);
  double responseSurfEval = responseSurf.eval(v);
  std::cout << "f(v)  = " << f(v) << std::endl;
  std::cout << "rI(v) = " << responseSurfEval << std::endl;

  std::cout << "err: " << responseSurf.l2Error(objectiveFunc) << std::endl;

  return 0;
}
