// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/optimization/activeSubspaces/ASMatrixGradientMC.hpp>
#include <sgpp/optimization/activeSubspaces/ASMatrixNakBspline.hpp>
#include <sgpp/optimization/activeSubspaces/ASResponseSurfaceNakBspline.hpp>
#include <sgpp/optimization/activeSubspaces/SparseGridResponseSurfaceNakBspline.hpp>
#include <sgpp/optimization/function/scalar/WrapperScalarFunction.hpp>
#include <sgpp/optimization/tools/Printer.hpp>

/*
 * Code to compare different methods of
 * 	finding active subspaces
 * 		known gradient and Monte Carlo
 * 		[differenc quotient and Monte Carlo]
 * 		regular B-spline interpolant and its exact Gauss integral
 * 		adaptive B-spline interpolant and its exact Gauss integral
 * 	and building the response surface
 * 		regression from the detection points
 * 		[by averaging over new adaptive grid]
 *
 * 	[...] stuff not yet available
 *
 */

// objective function

class objectiveFunction {
 public:
  objectiveFunction() {}
  ~objectiveFunction() {}
  static double f(sgpp::base::DataVector v) { return v[0] * v[0] + v[0] * v[1]; }
  static double df(sgpp::base::DataVector v, sgpp::base::DataVector& gradient) {
    gradient.resizeZero(2);
    gradient[0] = 2 * v[0] + v[1];
    gradient[1] = v[0];
    return v[0] * v[0] + v[0] * v[1];
  }
  sgpp::optimization::WrapperScalarFunction getObjectiveFunction() {
    size_t numDim = 2;
    return sgpp::optimization::WrapperScalarFunction(numDim, f);
  }
  sgpp::optimization::WrapperScalarFunctionGradient getObjectiveFunctionGradient() {
    size_t numDim = 2;
    return sgpp::optimization::WrapperScalarFunctionGradient(numDim, df);
  }
  Eigen::MatrixXd getActiveSubspaceMatrix() {
    Eigen::MatrixXd C(2, 2);
    C(0, 0) = 2.666666666666666;
    C(0, 1) = 0.916666666666666;
    C(1, 0) = 0.916666666666666;
    C(1, 1) = 0.333333333333333;
    return C;
  }
  Eigen::VectorXd getEigenValues() {
    Eigen::VectorXd e(2);
    e(0) = 0.016292182102929;
    e(1) = 2.983707817897071;
    return e;
  }
  Eigen::MatrixXd getEigenVectors() {
    Eigen::MatrixXd v(2, 2);
    v(0, 0) = -0.326865156584162;
    v(0, 1) = 0.945070986440284;
    v(1, 0) = 0.945070986440284;
    v(1, 1) = 0.326865156584162;
    return v;
  }
};

class objectiveFunctionExp {
 public:
  objectiveFunctionExp() {}
  ~objectiveFunctionExp() {}
  static double f(sgpp::base::DataVector v) { return exp(0.7 * v[0] + 0.3 * v[1]); }
  static double df(sgpp::base::DataVector v, sgpp::base::DataVector& gradient) {
    gradient.resizeZero(2);
    gradient[0] = 0.7 * exp(0.7 * v[0] + 0.3 * v[1]);
    gradient[1] = 0.3 * exp(0.7 * v[0] + 0.3 * v[1]);
    return exp(0.7 * v[0] + 0.3 * v[1]);
  }
  sgpp::optimization::WrapperScalarFunction getObjectiveFunction() {
    size_t numDim = 2;
    return sgpp::optimization::WrapperScalarFunction(numDim, f);
  }
  sgpp::optimization::WrapperScalarFunctionGradient getObjectiveFunctionGradient() {
    size_t numDim = 2;
    return sgpp::optimization::WrapperScalarFunctionGradient(numDim, df);
  }
  Eigen::MatrixXd getActiveSubspaceMatrix() {
    Eigen::MatrixXd C(2, 2);
    C(0, 0) = 1.46518;
    C(0, 1) = 0.627934;
    C(1, 0) = 0.627934;
    C(1, 1) = 0.269115;
    return C;
  }
  Eigen::VectorXd getEigenValues() {
    Eigen::VectorXd e(2);
    e(0) = 0.000000465517237;
    e(1) = 1.734294534482763;
    return e;
  }
  Eigen::MatrixXd getEigenVectors() {
    Eigen::MatrixXd v(2, 2);
    v(0, 0) = -0.393919252891044;
    v(0, 1) = 0.919145049598681;
    v(1, 0) = 0.919145049598681;
    v(1, 1) = 0.393919252891044;
    return v;
  }
};

class objectiveFunctionPoly {
 public:
  objectiveFunctionPoly() {}
  ~objectiveFunctionPoly() {}
  static double f(sgpp::base::DataVector v) { return sin(10 * v[0] + 1 * v[1] + 0 * v[2]); }
  static double df(sgpp::base::DataVector v, sgpp::base::DataVector& gradient) {
    gradient.resizeZero(3);
    gradient[0] = cos(10 * v[0] + 1 * v[1] + 0 * v[2]) * 10;
    gradient[1] = cos(10 * v[0] + 1 * v[1] + 0 * v[2]) * 1;
    gradient[2] = cos(10 * v[0] + 1 * v[1] + 0 * v[2]) * 0;
    return sin(10 * v[0] + 1 * v[1] + 0 * v[2]);
  }
  sgpp::optimization::WrapperScalarFunction getObjectiveFunction() {
    size_t numDim = 3;
    return sgpp::optimization::WrapperScalarFunction(numDim, f);
  }
  sgpp::optimization::WrapperScalarFunctionGradient getObjectiveFunctionGradient() {
    size_t numDim = 3;
    return sgpp::optimization::WrapperScalarFunctionGradient(numDim, df);
  }
};

class objectiveFunctionJeff {
 public:
  objectiveFunctionJeff() {}
  ~objectiveFunctionJeff() {}
  static double f(sgpp::base::DataVector v) {
    double p = 2;
    size_t n = 1;
    double partsum = 0.0;
    for (size_t j = 0; j < n - 1; j++) {
      partsum += std::pow(v[j], p - 1);
    }
    return std::pow(v.sum(), p) + partsum;
  }
  sgpp::optimization::WrapperScalarFunction getObjectiveFunction() {
    size_t numDim = 6;
    return sgpp::optimization::WrapperScalarFunction(numDim, f);
  }
};

class objectiveFunctionMaximumAS {
 public:
  objectiveFunctionMaximumAS() {}
  ~objectiveFunctionMaximumAS() {}
  static double f(sgpp::base::DataVector v) { return sin(10 * v[0] + v[1]); }
  sgpp::optimization::WrapperScalarFunction getObjectiveFunction() {
    size_t numDim = 2;
    return sgpp::optimization::WrapperScalarFunction(numDim, f);
  }
};

int main() {
  sgpp::optimization::Printer::getInstance().setVerbosity(-1);
  size_t degree = 3;
  // active subspace specifier
  size_t n = 1;
  size_t numMCErrorPoints = 1000;
  objectiveFunctionMaximumAS objectiveFuncInstance;
  sgpp::optimization::WrapperScalarFunction objectiveFunc =
      objectiveFuncInstance.getObjectiveFunction();
  //  sgpp::optimization::WrapperScalarFunctionGradient objectiveFuncGradient =
  //      objectiveFuncInstance.getObjectiveFunctionGradient();
  sgpp::base::GridType gridType = sgpp::base::GridType::NakBsplineBoundary;

  size_t numDim = objectiveFunc.getNumberOfParameters();
  for (size_t level = 1; level < 5; level++) {
    // regular sparse grid interpolant
    auto regularResponseSurf =
        std::make_shared<sgpp::optimization::SparseGridResponseSurfaceNakBspline>(
            numDim, objectiveFunc, gridType, degree);
    regularResponseSurf->createRegularResponseSurface(level);
    size_t numberOfPointsAccordingToLevel = regularResponseSurf->getSize();
    std::cout << "num points: " << numberOfPointsAccordingToLevel << std::endl;
    std::cout << "Reg  | l2 err: " << regularResponseSurf->l2Error(objectiveFunc, numMCErrorPoints)
              << std::endl;

    // surplus adaptive sparse grid interpolant
    auto adaptiveResponseSurf =
        std::make_shared<sgpp::optimization::SparseGridResponseSurfaceNakBspline>(
            numDim, objectiveFunc, gridType, degree);
    adaptiveResponseSurf->createSurplusAdaptiveResponseSurface(numberOfPointsAccordingToLevel, 1);
    std::cout << "Ad   | l2 err: " << adaptiveResponseSurf->l2Error(objectiveFunc, numMCErrorPoints)
              << std::endl;

    // adaptive B-spline matrix - regular Bspline response surface from detection points
    sgpp::optimization::ASMatrixNakBspline ASM_Bspline(objectiveFunc, gridType, degree);
    ASM_Bspline.buildAdaptiveInterpolant(numberOfPointsAccordingToLevel);
    ASM_Bspline.createMatrixGauss();
    ASM_Bspline.evDecompositionForSymmetricMatrices();
    Eigen::VectorXd eigenvalues_adaptiveBspline = ASM_Bspline.getEigenvalues();
    Eigen::MatrixXd eigenvectors_adaptiveBspline = ASM_Bspline.getEigenvectors();
    Eigen::MatrixXd W1_adaptiveBspline = ASM_Bspline.getTransformationMatrix(n);
    Eigen::MatrixXd C_adaptiveBspline = ASM_Bspline.getMatrix();
    sgpp::base::DataMatrix evaluationPoints_adaptiveBspline = ASM_Bspline.getEvaluationPoints();
    sgpp::base::DataVector functionValues_adaptiveBspline = ASM_Bspline.getFunctionValues();
    auto responseSurf_regularBspline =
        std::make_shared<sgpp::optimization::ASResponseSurfaceNakBspline>(
            numDim, W1_adaptiveBspline, gridType, degree);
    responseSurf_regularBspline->createRegularReducedSurfaceFromDetectionPoints(
        evaluationPoints_adaptiveBspline, functionValues_adaptiveBspline, level);
    std::cout << "B    | l2 err: "
              << responseSurf_regularBspline->l2Error(objectiveFunc, numMCErrorPoints)
              //              << " matrix err: "
              //              << (C_adaptiveBspline -
              //              objectiveFuncInstance.getActiveSubspaceMatrix()).norm()
              //              << " eigenval err: "
              //              << (eigenvalues_adaptiveBspline -
              //              objectiveFuncInstance.getEigenValues()).norm()
              //              << " eigenvec err: "
              //              << (eigenvectors_adaptiveBspline -
              //              objectiveFuncInstance.getEigenVectors()).norm()
              << std::endl;

    // adaptive B-spline matrix - regular Bspline response surface
    size_t maxNumGridPoints = functionValues_adaptiveBspline.getSize();
    size_t initialLevel = 1;
    auto responseSurf_regularPinvBspline =
        std::make_shared<sgpp::optimization::ASResponseSurfaceNakBspline>(
            numDim, W1_adaptiveBspline, gridType, degree);
    responseSurf_regularPinvBspline->createAdaptiveReducedSurfaceWithPseudoInverse(
        maxNumGridPoints, objectiveFunc, initialLevel);
    std::cout << "B r  | l2 err: "
              << responseSurf_regularPinvBspline->l2Error(objectiveFunc, numMCErrorPoints)
              << std::endl;

    // adaptive B-spline matrix - regular Bspline response surface
    auto responseSurf_adaptivePinvBspline =
        std::make_shared<sgpp::optimization::ASResponseSurfaceNakBspline>(
            numDim, W1_adaptiveBspline, gridType, degree);
    responseSurf_adaptivePinvBspline->createAdaptiveReducedSurfaceWithPseudoInverse(
        maxNumGridPoints, objectiveFunc, initialLevel);
    std::cout << "B a  | l2 err: "
              << responseSurf_adaptivePinvBspline->l2Error(objectiveFunc, numMCErrorPoints)
              << std::endl;

    // MC matrix with finite differences - regular Bspline response surface
    sgpp::optimization::ASMatrixGradientMC ASM_GradientMCFD(objectiveFunc);
    ASM_GradientMCFD.createMatrixMonteCarloFiniteDifference(numberOfPointsAccordingToLevel);
    ASM_GradientMCFD.evDecompositionForSymmetricMatrices();
    Eigen::VectorXd eigenvalues_gradientMCFD = ASM_GradientMCFD.getEigenvalues();
    Eigen::MatrixXd eigenvectors_gradientMCFD = ASM_GradientMCFD.getEigenvectors();
    Eigen::MatrixXd W1_gradientMCFD = ASM_GradientMCFD.getTransformationMatrix(n);
    Eigen::MatrixXd C_gradientMCFD = ASM_GradientMCFD.getMatrix();
    sgpp::base::DataMatrix evaluationPoints_gradientMCFD = ASM_GradientMCFD.getEvaluationPoints();
    sgpp::base::DataVector functionValues_gradientMCFD = ASM_GradientMCFD.getFunctionValues();
    auto responseSurf_gradientMCFD =
        std::make_shared<sgpp::optimization::ASResponseSurfaceNakBspline>(numDim, W1_gradientMCFD,
                                                                          gridType, degree);
    responseSurf_gradientMCFD->createRegularReducedSurfaceFromDetectionPoints(
        evaluationPoints_gradientMCFD, functionValues_gradientMCFD, level);
    std::cout
        << "MC FD| l2 err: "
        << responseSurf_gradientMCFD->l2Error(objectiveFunc, numMCErrorPoints)
        //              << " matrix err: "
        //              << (C_gradientMC - objectiveFuncInstance.getActiveSubspaceMatrix()).norm()
        //              << " eigenval err: "
        //              << (eigenvalues_gradientMC - objectiveFuncInstance.getEigenValues()).norm()
        //              << " eigenvec err: "
        //              << (eigenvectors_gradientMC -
        //              objectiveFuncInstance.getEigenVectors()).norm()
        << std::endl;

    // MC matrix with gradient - regular Bspline response surface
    //    sgpp::optimization::ASMatrixGradientMC ASM_GradientMC(objectiveFunc);
    //    ASM_GradientMC.createMatrixMonteCarlo(numberOfPointsAccordingToLevel,
    //    objectiveFuncGradient); ASM_GradientMC.evDecompositionForSymmetricMatrices();
    //    Eigen::VectorXd eigenvalues_gradientMC = ASM_GradientMC.getEigenvalues();
    //    Eigen::MatrixXd eigenvectors_gradientMC = ASM_GradientMC.getEigenvectors();
    //    Eigen::MatrixXd W1_gradientMC = ASM_GradientMC.getTransformationMatrix(n);
    //    Eigen::MatrixXd C_gradientMC = ASM_GradientMC.getMatrix();
    //    sgpp::base::DataMatrix evaluationPoints_gradientMC = ASM_GradientMC.getEvaluationPoints();
    //    sgpp::base::DataVector functionValues_gradientMC = ASM_GradientMC.getFunctionValues();
    //    auto responseSurf_gradientMC =
    //        std::make_shared<sgpp::optimization::ASResponseSurfaceNakBspline>(W1_gradientMC,
    //        gridType,
    //                                                                          degree);
    //    responseSurf_gradientMC->createRegularReducedSurfaceFromDetectionPoints(
    //        evaluationPoints_gradientMC, functionValues_gradientMC, level);
    //    std::cout
    //        << "MC   | l2 err: "
    //        << responseSurf_gradientMC->l2Error(objectiveFunc, numMCErrorPoints)
    //        //              << " matrix err: "
    //        //              << (C_gradientMC -
    //        objectiveFuncInstance.getActiveSubspaceMatrix()).norm()
    //        //              << " eigenval err: "
    //        //              << (eigenvalues_gradientMC -
    //        objectiveFuncInstance.getEigenValues()).norm()
    //        //              << " eigenvec err: "
    //        //              << (eigenvectors_gradientMC -
    //        //              objectiveFuncInstance.getEigenVectors()).norm()
    //        << std::endl;

    std::cout << "EVal:\n" << eigenvalues_adaptiveBspline << std::endl;
  }
  return 0;
}
