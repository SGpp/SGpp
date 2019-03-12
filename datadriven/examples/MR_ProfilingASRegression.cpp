// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/tools/SGppStopwatch.hpp>
#include <sgpp/datadriven/activeSubspaces/ASMatrixBsplineData.hpp>
#include <sgpp/datadriven/activeSubspaces/ASResponseSurfaceNakBspline.hpp>
#include <sgpp/optimization/function/scalar/WrapperScalarFunction.hpp>
#include <sgpp/optimization/tools/Printer.hpp>
#include <sgpp/optimization/tools/RandomNumberGenerator.hpp>

class objectiveFunctionVarious {
 public:
  objectiveFunctionVarious() {}
  ~objectiveFunctionVarious() {}
  static double f(sgpp::base::DataVector v) {
    return exp(0.7 * v[0] + 0.3 * v[1]);
    //    return exp(0.01 * v[0] + 0.02 * v[1] + 0.03 * v[2] + 0.04 * v[3] + 0.05 * v[4] +
    //    0.06 * v[5] +  0.07 * v[6] + 0.08 * v[7]);
  }
  sgpp::optimization::WrapperScalarFunction getObjectiveFunction() {
    size_t numDim = 2;
    return sgpp::optimization::WrapperScalarFunction(numDim, f);
  }
};

// Determine the active subspaces of a function given only by data. This is used for profiling.
int main() {
  sgpp::base::SGppStopwatch mainWatch;
  sgpp::base::SGppStopwatch watch;
  mainWatch.start();
  watch.start();
  sgpp::optimization::Printer::getInstance().setVerbosity(-1);
  size_t numDim = 2;
  size_t degree = 3;
  size_t n = 1;  // active subspace specifier
  size_t numInterpolPoints = 300;
  objectiveFunctionVarious objectiveFuncInstance;
  auto objectiveFunc = std::make_shared<sgpp::optimization::WrapperScalarFunction>(
      objectiveFuncInstance.getObjectiveFunction());
  sgpp::base::GridType gridType = sgpp::base::GridType::NakBsplineExtended;

  size_t numDataPoints = 10000;
  sgpp::base::DataMatrix evaluationPoints(numDim, numDataPoints);
  sgpp::base::DataVector functionValues(numDataPoints);
  sgpp::base::DataVector point(numDim);
  for (size_t i = 0; i < numDataPoints; i++) {
    sgpp::optimization::RandomNumberGenerator::getInstance().getUniformRV(point, 0.0, 1.0);
    evaluationPoints.setColumn(i, point);
    functionValues[i] = objectiveFunc->eval(point);
  }

  sgpp::datadriven::ASMatrixBsplineData ASM_Bspline(evaluationPoints, functionValues, gridType,
                                                    degree);
  ASM_Bspline.buildAdaptiveInterpolant(numInterpolPoints);
  std::cout << "constructing interpolant: " << watch.stop() << "s\n";
  watch.start();
  ASM_Bspline.createMatrixGauss();
  ASM_Bspline.evDecompositionForSymmetricMatrices();

  Eigen::VectorXd eigenvalues_adaptiveBspline = ASM_Bspline.getEigenvalues();
  Eigen::MatrixXd eigenvectors_adaptiveBspline = ASM_Bspline.getEigenvectors();
  Eigen::MatrixXd W1 = ASM_Bspline.getTransformationMatrix(n);
  Eigen::MatrixXd C_adaptiveBspline = ASM_Bspline.getMatrix();
  sgpp::base::DataMatrix evaluationPoints_adaptiveBspline = ASM_Bspline.getEvaluationPoints();
  sgpp::base::DataVector functionValues_adaptiveBspline = ASM_Bspline.getFunctionValues();

  std::cout << "EV decomposition: " << watch.stop() << "s\n";
  watch.start();
  auto responseSurf =
      std::make_shared<sgpp::datadriven::ASResponseSurfaceNakBspline>(W1, gridType, degree);
  responseSurf->createAdaptiveReducedSurfaceFromData(numInterpolPoints, evaluationPoints,
                                                     functionValues, 1, 3);
  std::cout << "interpolation error: " << responseSurf->l2Error(objectiveFunc, 1000) << "\n";
  std::cout << "constructing reduced interpolant: " << watch.stop() << "s\n";
  watch.start();

  //  double integral = responseSurf->getSplineBasedIntegral();
  double integral = responseSurf->getApproximateSplineBasedIntegral();
  std::cout << "integral = " << integral << "\n";
  std::cout << "integration: " << watch.stop() << "s\n";
  std::cout << "total runtime: " << mainWatch.stop() << "\n";
  return 0;
}
