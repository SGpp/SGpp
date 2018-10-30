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
  //    return exp(0.7 * v[0] + 0.3 * v[1]);
  //  return sin(100 * v[0] + 10 * v[1] + 1 * v[2] + 0 * v[3]);
  return exp(5 * v[0] + v[1]);
}
double df(sgpp::base::DataVector v, sgpp::base::DataVector& gradient) {
  gradient.resizeZero(4);
  gradient[0] = cos(100 * v[0] + 1 * v[1]) * 100;
  gradient[1] = cos(100 * v[0] + 1 * v[1]) * 10;
  gradient[2] = cos(100 * v[0] + 1 * v[1]) * 1;
  gradient[3] = cos(100 * v[0] + 1 * v[1]) * 0;
  return sin(100 * v[0] + 10 * v[1] + 1 * v[2] + 0 * v[3]);
}

int main() {
  sgpp::optimization::Printer::getInstance().setVerbosity(-1);
  size_t numDim = 2;
  size_t degree = 3;
  // active subspace specifier
  size_t n = 1;

  auto objectiveFunc = std::make_shared<sgpp::optimization::WrapperScalarFunction>(numDim, f);
  auto objectiveFuncGradient =
      std::make_shared<sgpp::optimization::WrapperScalarFunctionGradient>(numDim, df);
  sgpp::base::GridType gridType = sgpp::base::GridType::NakBsplineBoundary;
  sgpp::optimization::ASMatrixNakBspline ASM(objectiveFunc, gridType, degree);
  size_t maxNumGridPointsMatrix = 200;
  size_t initialLevel = 1;
  ASM.buildAdaptiveInterpolant(maxNumGridPointsMatrix, initialLevel);
  ASM.createMatrixGauss();
  ASM.evDecompositionForSymmetricMatrices();

  Eigen::VectorXd eigenvalues = ASM.getEigenvalues();
  std::cout << "EVal:\n" << eigenvalues << std::endl;
  Eigen::MatrixXd eigenvectors = ASM.getEigenvectors();
  std::cout << "EVec:\n" << eigenvectors << std::endl;
  Eigen::MatrixXd W1 = ASM.getTransformationMatrix(n);
  std::cout << "W1\n" << W1 << std::endl;
  std::cout << "-----------------------------" << std::endl;

  sgpp::base::DataMatrix evaluationPoints = ASM.getEvaluationPoints();
  sgpp::base::DataVector functionValues = ASM.getFunctionValues();

  for (size_t responseLevel = 1; responseLevel < 6; responseLevel++) {
    sgpp::optimization::ASResponseSurfaceNakBspline responseSurf(numDim, W1, gridType, degree);
    //  size_t maxNumGridPointsResponseSurface = 500;
    //  responseSurf.createAdaptiveReducedSurfaceWithPseudoInverse(maxNumGridPointsResponseSurface,
    //                                                             objectiveFunc, initialLevel);
    responseSurf.createRegularReducedSurfaceWithPseudoInverse(responseLevel, objectiveFunc);

    std::cout << responseLevel << " err: " << responseSurf.l2Error(objectiveFunc) << std::endl;
    //    sgpp::base::DataVector v(numDim, 1);
    //    double responseSurfEval = responseSurf.eval(v);
    //    std::cout << "f(v)  = " << f(v) << std::endl;
    //    std::cout << "rI(v) = " << responseSurfEval << std::endl;
  }

  //  sgpp::base::DataVector v(numDim, 0.3371);
  //  double responseSurfEval = responseSurf.eval(v);
  //  std::cout << "f(v)  = " << f(v) << std::endl;
  //  std::cout << "rI(v) = " << responseSurfEval << std::endl;

  //  std::cout << "--------- MC ---------" << std::endl;
  //  size_t numMCpoints = 100;
  //  sgpp::optimization::ASMatrixGradientMC ASM_GradientMC(objectiveFunc);
  //  ASM_GradientMC.createMatrixMonteCarlo(numMCpoints, objectiveFuncGradient);
  //  ASM_GradientMC.evDecompositionForSymmetricMatrices();
  //  Eigen::VectorXd eigenvalues_gradientMC = ASM_GradientMC.getEigenvalues();
  //  Eigen::MatrixXd eigenvectors_gradientMC = ASM_GradientMC.getEigenvectors();
  //  Eigen::MatrixXd W1_gradientMC = ASM_GradientMC.getTransformationMatrix(n);
  //  std::cout << "EVal\n" << eigenvalues_gradientMC << std::endl;
  //  std::cout << "EVec\n" << eigenvectors_gradientMC << std::endl;
  //  std::cout << "W1\n" << W1_gradientMC << std::endl;
  //  sgpp::base::DataMatrix evaluationPoints_gradientMC = ASM_GradientMC.getEvaluationPoints();
  //  sgpp::base::DataVector functionValues_gradientMC = ASM_GradientMC.getFunctionValues();
  //  auto responseSurf_gradientMC =
  //  std::make_shared<sgpp::optimization::ASResponseSurfaceNakBspline>(dim,
  //      W1_gradientMC, gridType, degree);
  //  responseSurf_gradientMC->createAdaptiveReducedSurfaceWithPseudoInverse(
  //      maxNumGridPointsResponseSurface, objectiveFunc, initialLevel);
  //  std::cout << "err: " << responseSurf_gradientMC->l2Error(objectiveFunc) << std::endl;
  return 0;
}
