// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/tools/SGppStopwatch.hpp>
#include <sgpp/datadriven/activeSubspaces/ASMatrixBsplineAnalytic.hpp>
#include <sgpp/datadriven/activeSubspaces/ASResponseSurfaceNakBspline.hpp>
#include <sgpp/optimization/function/scalar/WrapperScalarFunction.hpp>
#include <sgpp/optimization/tools/Printer.hpp>

class objectiveFunctionVarious {
 public:
  objectiveFunctionVarious() {}
  ~objectiveFunctionVarious() {}
  static double f(sgpp::base::DataVector v) {
    return exp(0.01 * v[0] + 0.02 * v[1] + 0.03 * v[2] + 0.04 * v[3] + 0.05 * v[4] + 00.06 * v[5]);
  }
  sgpp::optimization::WrapperScalarFunction getObjectiveFunction() {
    size_t numDim = 6;
    return sgpp::optimization::WrapperScalarFunction(numDim, f);
  }
};

int main() {
  sgpp::base::SGppStopwatch mainWatch;
  sgpp::base::SGppStopwatch watch;
  mainWatch.start();
  watch.start();
  sgpp::optimization::Printer::getInstance().setVerbosity(-1);
  size_t degree = 3;
  size_t n = 1;  // active subspace specifier
  size_t numInterpolPoints = 1000;
  objectiveFunctionVarious objectiveFuncInstance;
  auto objectiveFunc = std::make_shared<sgpp::optimization::WrapperScalarFunction>(
      objectiveFuncInstance.getObjectiveFunction());
  sgpp::base::GridType gridType = sgpp::base::GridType::NakBsplineModified;

  sgpp::datadriven::ASMatrixBsplineAnalytic ASM(objectiveFunc, gridType, degree);
  ASM.buildAdaptiveInterpolant(numInterpolPoints);
  std::cout << "constructing interpolant: " << watch.stop() << "s\n";
  watch.start();
  ASM.createMatrixGauss();
  ASM.evDecompositionForSymmetricMatrices();

  Eigen::VectorXd eigenvalues_adaptiveBspline = ASM.getEigenvalues();
  Eigen::MatrixXd eigenvectors_adaptiveBspline = ASM.getEigenvectors();
  Eigen::MatrixXd W1 = ASM.getTransformationMatrix(n);
  Eigen::MatrixXd C_adaptiveBspline = ASM.getMatrix();
  sgpp::base::DataMatrix evaluationPoints_adaptiveBspline = ASM.getEvaluationPoints();
  sgpp::base::DataVector functionValues_adaptiveBspline = ASM.getFunctionValues();

  std::cout << "EV decomposition: " << watch.stop() << "s\n";
  watch.start();
  //  auto responseSurf =
  //      std::make_shared<sgpp::datadriven::ASResponseSurfaceNakBspline>(W1, gridType, degree);

  sgpp::base::GridType responseGridType = sgpp::base::GridType::NakBsplineBoundary;
  auto responseSurf = ASM.getResponseSurfaceInstance(n, responseGridType, degree);
  responseSurf.createAdaptiveReducedSurfaceWithPseudoInverse(numInterpolPoints, objectiveFunc, 1,
                                                             3);

  std::cout << "constructing reduced interpolant: " << watch.stop() << "s\n";
  watch.start();

  //  double integral = responseSurf->getSplineBasedIntegral();
  size_t quadOrder = 7;
  double integral = responseSurf.getSplineBasedIntegral(quadOrder);
  std::cout << "integral = " << integral << "\n";
  std::cout << "integration: " << watch.stop() << "s\n";
  std::cout << "total runtime: " << mainWatch.stop() << "\n";
  return 0;
}
