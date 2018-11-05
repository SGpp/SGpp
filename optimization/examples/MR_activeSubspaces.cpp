// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/tools/SGppStopwatch.hpp>
#include <sgpp/optimization/activeSubspaces/ASMatrixGradientMC.hpp>
#include <sgpp/optimization/activeSubspaces/ASMatrixNakBspline.hpp>
#include <sgpp/optimization/activeSubspaces/ASResponseSurfaceNakBspline.hpp>
#include <sgpp/optimization/function/scalar/WrapperScalarFunction.hpp>
#include <sgpp/optimization/tools/Printer.hpp>

#include <cmath>
#include <sstream>

void initialize(int argc, char* argv[], bool& adaptive, bool& print, size_t& degree, size_t& level,
                size_t& numPoints, sgpp::base::GridType& gridType) {
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if ((arg == "-r") || (arg == "--regular")) {
      adaptive = 0;
    } else if ((arg == "-a") || (arg == "--adaptive")) {
      adaptive = 1;
    } else if ((arg == "-p") || (arg == "--print")) {
      print = 1;
    } else if ((arg == "-d") || (arg == "--degree")) {
      std::string degreeString = argv[++i];
      std::stringstream degreeStream(degreeString);
      degreeStream >> degree;
    } else if ((arg == "-l") || (arg == "--level")) {
      std::string levelString = argv[++i];
      std::stringstream levelStream(levelString);
      levelStream >> level;
    } else if ((arg == "-n") || (arg == "--numPoints")) {
      std::string numPointsString = argv[++i];
      std::stringstream numPointsStream(numPointsString);
      numPointsStream >> numPoints;
    } else if ((arg == "-g") || (arg == "--gridType")) {
      std::string gridTypeString = argv[++i];
      gridType = sgpp::base::Grid::stringToGridType(gridTypeString);
    }
  }
}

// function from basic.ipynb in active_subspaces/tutorials
// input: v in  [0,1]^10 must be transformed to actual domain given through lowerBounds, upperBounds
double wing(sgpp::base::DataVector v) {
  std::vector<double> lowerBounds{150, 220, 6, -10, 16, .5, .08, 2.5, 1700, .025};
  std::vector<double> upperBounds{200, 300, 10, 10, 45, 1, .18, 6, 2500, .08};
  double Sw = v[0] * (upperBounds[0] - lowerBounds[0]) + lowerBounds[0];
  double Wfw = v[1] * (upperBounds[1] - lowerBounds[1]) + lowerBounds[1];
  double A = v[2] * (upperBounds[2] - lowerBounds[2]) + lowerBounds[2];
  double Lambda = v[3] * (upperBounds[3] - lowerBounds[3]) + lowerBounds[3];
  double q = v[4] * (upperBounds[4] - lowerBounds[4]) + lowerBounds[4];
  double lambda = v[5] * (upperBounds[5] - lowerBounds[5]) + lowerBounds[5];
  double tc = v[6] * (upperBounds[6] - lowerBounds[6]) + lowerBounds[6];
  double Nz = v[7] * (upperBounds[7] - lowerBounds[7]) + lowerBounds[7];
  double Wdg = v[8] * (upperBounds[8] - lowerBounds[8]) + lowerBounds[8];
  double Wp = v[9] * (upperBounds[9] - lowerBounds[9]) + lowerBounds[9];

  double res1 = 0.036 * std::pow(Sw, 0.758);
  double res2 = std::pow(Wfw, 0.0035);
  double res3 = std::pow(A, 0.6);
  double res4 = std::pow(q, 0.006);
  double res5 = std::pow(lambda, 0.04);
  double res6 = std::pow(100 * tc, -0.3);
  double res7 = std::pow(Nz * Wdg, 0.49);
  double res8 = std::pow(std::cos(Lambda), -0.9);
  double res9 = Sw * Wp;

  // res8 can equal nan. Set to 1 if that happens.
  // This is not physically motivated and probably leads to wrong results!!!
  if (std::isnan(res8)) {
    res8 = 1;
  }
  double res = res1 * res2 * res3 * res4 * res5 * res6 * res7 * res8 + res9;
  return res;
}

double f(sgpp::base::DataVector v) {
  //  return exp(0.7 * v[0] + 0.3 * v[1]);
  //  return sin(100 * v[0] + 10 * v[1] + 1 * v[2] + 0 * v[3]);
  //  return exp(5 * v[0] + v[1]);
  //  return sin(v[0] + v[1]) * cos(v[0] + v[1]);
  return wing(v);
}
double df(sgpp::base::DataVector v, sgpp::base::DataVector& gradient) {
  gradient.resizeZero(4);
  gradient[0] = cos(100 * v[0] + 1 * v[1]) * 100;
  gradient[1] = cos(100 * v[0] + 1 * v[1]) * 10;
  gradient[2] = cos(100 * v[0] + 1 * v[1]) * 1;
  gradient[3] = cos(100 * v[0] + 1 * v[1]) * 0;
  return sin(100 * v[0] + 10 * v[1] + 1 * v[2] + 0 * v[3]);
}

int main(int argc, char* argv[]) {
  // command line parameters
  bool adaptive = 0;
  bool print = 0;
  size_t degree = 3;
  size_t maxLevel = 9;
  size_t maxNumPoints = 100;
  sgpp::base::GridType gridType = sgpp::base::GridType::NakBsplineModified;
  initialize(argc, argv, adaptive, print, degree, maxLevel, maxNumPoints, gridType);
  sgpp::optimization::Printer::getInstance().setVerbosity(-1);
  size_t numDim = 10;
  // active subspace specifier
  size_t n = 1;

  sgpp::base::SGppStopwatch watch;

  auto objectiveFunc = std::make_shared<sgpp::optimization::WrapperScalarFunction>(numDim, f);
  auto objectiveFuncGradient =
      std::make_shared<sgpp::optimization::WrapperScalarFunctionGradient>(numDim, df);
  sgpp::optimization::ASMatrixNakBspline ASM(objectiveFunc, gridType, degree);
  size_t maxNumGridPointsMatrix = 700;
  size_t initialLevel = 1;
  watch.start();
  ASM.buildAdaptiveInterpolant(maxNumGridPointsMatrix, initialLevel, 3);
  std::cout << watch.stop() << "s\n";
  watch.start();
  ASM.createMatrixGauss();
  std::cout << watch.stop() << "s\n";
  ASM.evDecompositionForSymmetricMatrices();

  Eigen::VectorXd eigenvalues = ASM.getEigenvalues();
  Eigen::MatrixXd eigenvectors = ASM.getEigenvectors();
  Eigen::MatrixXd W1 = ASM.getTransformationMatrix(n);
  if (print) {
    std::cout << "EVal:\n" << eigenvalues << std::endl;
    std::cout << "EVec:\n" << eigenvectors << std::endl;
    std::cout << "W1\n" << W1 << std::endl;
    std::cout << "-----------------------------" << std::endl;
  }

  sgpp::base::DataMatrix evaluationPoints = ASM.getEvaluationPoints();
  sgpp::base::DataVector functionValues = ASM.getFunctionValues();

  for (size_t responseLevel = 1; responseLevel <= maxLevel; responseLevel++) {
    sgpp::optimization::ASResponseSurfaceNakBspline responseSurf(numDim, W1, gridType, degree);
    //  size_t maxNumGridPointsResponseSurface = 500;
    //  responseSurf.createAdaptiveReducedSurfaceWithPseudoInverse(maxNumGridPointsResponseSurface,
    //                                                             objectiveFunc, initialLevel);
    responseSurf.createRegularReducedSurfaceWithPseudoInverse(responseLevel, objectiveFunc);

    std::cout << "l " << responseLevel << " " << responseSurf.getSize()
              << " pts err: " << responseSurf.l2Error(objectiveFunc, 3000) << "\n";
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
