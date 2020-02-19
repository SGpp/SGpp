// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/function/scalar/WrapperScalarFunction.hpp>
#include <sgpp/base/tools/Printer.hpp>
#include <sgpp/base/tools/SGppStopwatch.hpp>
#include <sgpp/datadriven/activeSubspaces/ASMatrixBsplineAnalytic.hpp>
#include <sgpp/datadriven/activeSubspaces/ASResponseSurfaceNakBspline.hpp>

class objectiveFunctionVarious {
 public:
  objectiveFunctionVarious() {}
  ~objectiveFunctionVarious() {}
  static double f(sgpp::base::DataVector v) { return sin(M_PI * v.sum()); }
  sgpp::base::WrapperScalarFunction getObjectiveFunction() {
    size_t numDim = 6;
    return sgpp::base::WrapperScalarFunction(numDim, f);
  }
};

// Calculate the integral of a function with a one-dim active subspace using the spline based
// method. This is used for profiling.
int main() {
  sgpp::base::SGppStopwatch mainWatch;
  sgpp::base::SGppStopwatch watch;
  mainWatch.start();
  watch.start();
  sgpp::base::Printer::getInstance().setVerbosity(-1);
  size_t degree = 3;
  size_t n = 1;  // active subspace specifier
  size_t numInterpolPoints = 2000;
  size_t numASM = numInterpolPoints;
  objectiveFunctionVarious objectiveFuncInstance;
  auto objectiveFunc = std::make_shared<sgpp::base::WrapperScalarFunction>(
      objectiveFuncInstance.getObjectiveFunction());
  sgpp::base::GridType gridType = sgpp::base::GridType::ModNakBspline;

  sgpp::datadriven::ASMatrixBsplineAnalytic ASM(objectiveFunc, gridType, degree);
  size_t initialLevel = 1;
  size_t numRefine = 3;  // 40;
  ASM.buildAdaptiveInterpolant(numASM, initialLevel, numRefine);
  std::cout << "constructing interpolant: " << watch.stop() << "s\n";
  watch.start();
  ASM.createMatrixGauss();
  ASM.evDecompositionForSymmetricMatrices();

  //  Eigen::VectorXd eigenvalues = ASM.getEigenvalues();
  //  Eigen::MatrixXd eigenvectors = ASM.getEigenvectors();
  //  Eigen::MatrixXd W1 = ASM.getTransformationMatrix(n);
  //  Eigen::MatrixXd C = ASM.getMatrix();
  //  sgpp::base::DataMatrix evaluationPoints = ASM.getEvaluationPoints();
  //  sgpp::base::DataVector functionValues = ASM.getFunctionValues();

  std::cout << "EV decomposition: " << watch.stop() << "s\n";
  watch.start();
  //  auto responseSurf =
  //      std::make_shared<sgpp::datadriven::ASResponseSurfaceNakBspline>(W1, gridType, degree);

  sgpp::base::GridType responseGridType = sgpp::base::GridType::NakBsplineBoundary;
  auto responseSurf = ASM.getResponseSurfaceInstance(n, responseGridType, degree);
  responseSurf.createAdaptiveReducedSurfaceWithPseudoInverse(numInterpolPoints, objectiveFunc,
                                                             initialLevel, numRefine);

  std::cout << "constructing reduced interpolant: " << watch.stop() << "s\n";
  std::cout << "l2 Error: " << responseSurf.l2Error(objectiveFunc, 10000) << "\n";
  watch.start();

  //  double integral = responseSurf->getSplineBasedIntegral();
  size_t quadOrder = 7;
  double integral = responseSurf.getSplineBasedIntegral(quadOrder);
  std::cout << "integration: " << watch.stop() << "s\n";
  std::cout << "integral = " << integral << "s\n";
  // only valid for sin5D
  //  std::cout << "integral error: " << fabs(integral - 32 / (M_PI * M_PI * M_PI * M_PI * M_PI))
  //            << "\n";
  std::cout << "total runtime: " << mainWatch.stop() << "\n";
  return 0;
}
