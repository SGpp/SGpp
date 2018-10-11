// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/optimization/activeSubspaces/ASMatrixGradientMC.hpp>
#include <sgpp/optimization/activeSubspaces/ASMatrixNakBspline.hpp>
#include <sgpp/optimization/activeSubspaces/ASResponseSurfaceNakBspline.hpp>
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

class objectiveFunction2 {
 public:
  objectiveFunction2() {}
  ~objectiveFunction2() {}
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

double l2InterpolationError(sgpp::optimization::WrapperScalarFunction func,
                            std::shared_ptr<sgpp::optimization::ASResponseSurface> reSurf,
                            size_t fullLevel) {
  double err = 0.0;
  size_t numDim = func.getNumberOfParameters();
  size_t dummyDegree = 3;
  std::shared_ptr<sgpp::base::Grid> grid;
  grid = std::make_shared<sgpp::base::NakBsplineBoundaryGrid>(numDim, dummyDegree);
  grid->getGenerator().full(fullLevel);
  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  for (size_t j = 0; j < gridStorage.getSize(); j++) {
    sgpp::base::GridPoint& gp = gridStorage.getPoint(j);
    sgpp::base::DataVector point(numDim);
    gp.getStandardCoordinates(point);
    err += std::pow(func.eval(point) - reSurf->eval(point), 2);
  }

  //  std::cout << "calculated error on " << grid->getSize() << " points." << std::endl;
  return sqrt(err);
}

int main() {
  sgpp::optimization::Printer::getInstance().setVerbosity(-1);
  size_t degree = 3;
  // active subspace specifier
  size_t n = 1;
  size_t responseLevel = 5;  // How to choose this? Adaptivity ?
  size_t errorGridLevel = 5;
  objectiveFunction2 objectiveFuncInstance;
  sgpp::optimization::WrapperScalarFunction objectiveFunc =
      objectiveFuncInstance.getObjectiveFunction();
  sgpp::optimization::WrapperScalarFunctionGradient objectiveFuncGradient =
      objectiveFuncInstance.getObjectiveFunctionGradient();
  sgpp::base::GridType gridType = sgpp::base::GridType::NakBsplineBoundary;

  std::vector<size_t> numPointsVector{10, 100, 300};
  for (auto& numPoints : numPointsVector) {
    std::cout << "num points: " << numPoints << std::endl;
    // adaptive B-spline
    sgpp::optimization::ASMatrixNakBspline ASM_Bspline(objectiveFunc, gridType, degree);
    ASM_Bspline.buildAdaptiveInterpolant(numPoints);
    ASM_Bspline.createMatrixGauss();
    ASM_Bspline.evDecompositionForSymmetricMatrices();
    Eigen::VectorXd eigenvalues_adaptiveBspline = ASM_Bspline.getEigenvalues();
    Eigen::MatrixXd eigenvectors_adaptiveBspline = ASM_Bspline.getEigenvectors();
    Eigen::MatrixXd W1_adaptiveBspline = ASM_Bspline.getTransformationMatrix(n);
    Eigen::MatrixXd C_adaptiveBspline = ASM_Bspline.getMatrix();
    sgpp::base::DataMatrix evaluationPoints_adaptiveBspline = ASM_Bspline.getEvaluationPoints();
    sgpp::base::DataVector functionValues_adaptiveBspline = ASM_Bspline.getFunctionValues();
    auto responseSurf_adaptiveBspline =
        std::make_shared<sgpp::optimization::ASResponseSurfaceNakBspline>(W1_adaptiveBspline,
                                                                          gridType, degree);
    responseSurf_adaptiveBspline->createRegularSurfaceFromDetectionPoints(
        evaluationPoints_adaptiveBspline, functionValues_adaptiveBspline, responseLevel);
    double error_adaptiveBspline =
        l2InterpolationError(objectiveFunc, responseSurf_adaptiveBspline, errorGridLevel);
    double matrixError_adaptiveBspline =
        (C_adaptiveBspline - objectiveFuncInstance.getActiveSubspaceMatrix()).norm();
    double eigenvalError_adaptiveBspline =
        (eigenvalues_adaptiveBspline - objectiveFuncInstance.getEigenValues()).norm();
    double eigenvecError_adaptiveBspline =
        (eigenvectors_adaptiveBspline - objectiveFuncInstance.getEigenVectors()).norm();
    std::cout << "B  | l2 err: " << error_adaptiveBspline
              << " matrix err: " << matrixError_adaptiveBspline
              << " eigenval err: " << eigenvalError_adaptiveBspline
              << " eigenvec err: " << eigenvecError_adaptiveBspline << std::endl;

    // MC
    sgpp::optimization::ASMatrixGradientMC ASM_GradientMC(objectiveFunc, objectiveFuncGradient);
    ASM_GradientMC.createMatrix(numPoints);
    ASM_GradientMC.evDecompositionForSymmetricMatrices();
    Eigen::VectorXd eigenvalues_gradientMC = ASM_GradientMC.getEigenvalues();
    Eigen::MatrixXd eigenvectors_gradientMC = ASM_GradientMC.getEigenvectors();
    Eigen::MatrixXd W1_gradientMC = ASM_GradientMC.getTransformationMatrix(n);
    Eigen::MatrixXd C_gradientMC = ASM_GradientMC.getMatrix();
    sgpp::base::DataMatrix evaluationPoints_gradientMC = ASM_GradientMC.getEvaluationPoints();
    sgpp::base::DataVector functionValues_gradientMC = ASM_GradientMC.getFunctionValues();
    auto responseSurf_gradientMC =
        std::make_shared<sgpp::optimization::ASResponseSurfaceNakBspline>(W1_gradientMC, gridType,
                                                                          degree);
    responseSurf_gradientMC->createRegularSurfaceFromDetectionPoints(
        evaluationPoints_gradientMC, functionValues_gradientMC, responseLevel);
    double error_gradientMC =
        l2InterpolationError(objectiveFunc, responseSurf_gradientMC, errorGridLevel);
    double matrixError_gradientMC =
        (C_gradientMC - objectiveFuncInstance.getActiveSubspaceMatrix()).norm();
    double eigenvalError_gradientMC =
        (eigenvalues_gradientMC - objectiveFuncInstance.getEigenValues()).norm();
    double eigenvecError_gradientMC =
        (eigenvectors_gradientMC - objectiveFuncInstance.getEigenVectors()).norm();
    std::cout << "MC | l2 err: " << error_gradientMC << " matrix err: " << matrixError_gradientMC
              << " eigenval err: " << eigenvalError_gradientMC
              << " eigenvec err: " << eigenvecError_gradientMC << std::endl;
  }
  return 0;
}
