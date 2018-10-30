// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp/base/grid/type/NakBsplineBoundaryGrid.hpp>
#include <sgpp/base/grid/type/NakBsplineModifiedGrid.hpp>
#include <sgpp/optimization/activeSubspaces/ASMatrixNakBspline.hpp>
#include <sgpp/optimization/activeSubspaces/ASResponseSurfaceNakBspline.hpp>
#include <sgpp/optimization/activeSubspaces/EigenFunctionalities.hpp>
#include <sgpp/optimization/activeSubspaces/GaussQuadrature.hpp>
#include <sgpp/optimization/function/scalar/InterpolantScalarFunction.hpp>
#include <sgpp/optimization/function/scalar/WrapperScalarFunction.hpp>
#include <sgpp/optimization/tools/Printer.hpp>

#include <functional>

double dummyFunction(sgpp::base::DataVector v) { return 777; }

double oneFunction(double x) { return 1; }
double xFunction(double x) { return x; }
double x6Function(double x) { return x * x * x * x * x * x; }

BOOST_AUTO_TEST_SUITE(testActiveSubspaces)

//#ifdef USE_EIGEN
BOOST_AUTO_TEST_CASE(testEigenFunctionalities) {
  sgpp::base::DataVector v(3);
  v[0] = 1;
  v[1] = -17.3;
  v[2] = 23;
  Eigen::VectorXd e = sgpp::optimization::DataVectorToEigen(v);
  //  std::cout << v.getSize() << " " << e.rows() << std::endl;
  //  std::cout << v[0] << " " << e(0) << std::endl;
  //  std::cout << v[1] << " " << e(1) << std::endl;
  //  std::cout << v[2] << " " << e(2) << std::endl;
  BOOST_CHECK_EQUAL(v.getSize(), e.rows());
  BOOST_CHECK_EQUAL(v[0], e(0));
  BOOST_CHECK_EQUAL(v[1], e(1));
  BOOST_CHECK_EQUAL(v[2], e(2));

  Eigen::VectorXd e2(3);
  e2(0) = -2;
  e2(1) = 11.1;
  e2(2) = 342;
  sgpp::base::DataVector v2 = sgpp::optimization::EigenToDataVector(e2);
  //  std::cout << v2.getSize() << " " << e2.rows() << std::endl;
  //  std::cout << v2[0] << " " << e2(0) << std::endl;
  //  std::cout << v2[1] << " " << e2(1) << std::endl;
  //  std::cout << v2[2] << " " << e2(2) << std::endl;
  BOOST_CHECK_EQUAL(v2.getSize(), e2.rows());
  BOOST_CHECK_EQUAL(v2[0], e2(0));
  BOOST_CHECK_EQUAL(v2[1], e2(1));
  BOOST_CHECK_EQUAL(v2[2], e2(2));
}
//#endif /* USE_EIGEN */

BOOST_AUTO_TEST_CASE(testQuad) {
  double epsilon = 1e-15;

  std::function<double(double)> func = oneFunction;
  size_t quadOrder = 1;
  double correctRes = 5;
  double quadResult = sgpp::optimization::gaussQuad(func, -3, 2, quadOrder);
  double err = fabs(correctRes - quadResult);
  //  std::cout << err << std::endl;
  BOOST_CHECK_SMALL(err, epsilon);

  func = xFunction;
  quadOrder = 1;
  correctRes = 50;
  quadResult = sgpp::optimization::gaussQuad(func, 0, 10, quadOrder);
  err = fabs(correctRes - quadResult);
  //  std::cout << err << std::endl;
  BOOST_CHECK_SMALL(err, epsilon);

  func = x6Function;
  quadOrder = 4;
  correctRes = 1.0 / 7.0;
  quadResult = sgpp::optimization::gaussQuad(func, 0, 1, quadOrder);
  err = fabs(correctRes - quadResult);
  //  std::cout << err << std::endl;
  BOOST_CHECK_SMALL(err, epsilon);
}

sgpp::base::DataVector interpolateRegular1D(
    std::shared_ptr<sgpp::optimization::WrapperScalarFunction> objectiveFunction,
    sgpp::base::GridType gridType, size_t degree, size_t level, sgpp::base::Grid* grid) {
  sgpp::optimization::Printer::getInstance().setVerbosity(-1);
  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  grid->getGenerator().regular(level);
  sgpp::base::DataVector functionValues(gridStorage.getSize());
  for (size_t i = 0; i < gridStorage.getSize(); i++) {
    sgpp::base::GridPoint& gp = gridStorage.getPoint(i);
    functionValues[i] =
        objectiveFunction->eval(sgpp::base::DataVector(1, gp.getStandardCoordinate(0)));
  }

  // solve linear system
  sgpp::base::DataVector alpha(functionValues.getSize());
  sgpp::optimization::HierarchisationSLE hierSLE(*grid);
  sgpp::optimization::sle_solver::Armadillo sleSolver;
  if (!sleSolver.solve(hierSLE, functionValues, alpha)) {
    std::cout << "TestASMatrixNakBsplineScalarProduct: Solving failed.\n";
  }
  return alpha;
}

double objectiveFunctionScalarProduct(sgpp::base::DataVector v) { return v[0] * v[0] * v[0]; }
double objectiveFunctionScalarProduct2(sgpp::base::DataVector v) { return sin(v[0]); }

BOOST_AUTO_TEST_CASE(testASMatrixNakBsplineBoundaryScalarProduct) {
  // interpolate function f = x^3  and calculate
  // \int f^2 dx
  // \int f'*f dx
  // \int f' * f' dx
  // with the scalar product routines of ASMatrixNakBspline

  size_t degree = 3;
  size_t level = 3;
  size_t numDim = 1;
  sgpp::base::GridType gridType = sgpp::base::GridType::NakBsplineBoundary;
  sgpp::base::Grid* grid;
  sgpp::base::Grid* grid2;
  if (gridType == sgpp::base::GridType::NakBspline) {
    grid = new sgpp::base::NakBsplineGrid(numDim, degree);
    grid2 = new sgpp::base::NakBsplineGrid(numDim, degree);
  } else if (gridType == sgpp::base::GridType::NakBsplineModified) {
    grid = new sgpp::base::NakBsplineModifiedGrid(numDim, degree);
    grid2 = new sgpp::base::NakBsplineModifiedGrid(numDim, degree);
  } else if (gridType == sgpp::base::GridType::NakBsplineBoundary) {
    grid = new sgpp::base::NakBsplineBoundaryGrid(numDim, degree);
    grid2 = new sgpp::base::NakBsplineBoundaryGrid(numDim, degree);
  } else {
    throw sgpp::base::generation_exception(
        "testASMatrixNakBsplineBoundaryScalarProducte: gridType not supported.");
  }

  // prepare gauss quadrature rule
  sgpp::base::DataVector coordinates, weights;
  sgpp::base::GaussLegendreQuadRule1D gauss;
  size_t quadOrder = static_cast<size_t>(std::ceil(static_cast<double>(degree) + 1.0 / 2.0));
  gauss.getLevelPointsAndWeightsNormalized(quadOrder, coordinates, weights);
  auto pCoordinates = std::make_shared<sgpp::base::DataVector>(coordinates);
  auto pWeights = std::make_shared<sgpp::base::DataVector>(weights);

  auto objectiveFunc = std::make_shared<sgpp::optimization::WrapperScalarFunction>(
      numDim, objectiveFunctionScalarProduct);
  auto objectiveFunc2 = std::make_shared<sgpp::optimization::WrapperScalarFunction>(
      numDim, objectiveFunctionScalarProduct2);
  sgpp::base::DataVector alpha = interpolateRegular1D(objectiveFunc, gridType, degree, level, grid);
  sgpp::base::DataVector alpha2 =
      interpolateRegular1D(objectiveFunc2, gridType, degree, level, grid2);
  sgpp::optimization::ASMatrixNakBspline ASM(objectiveFunc, gridType, degree);

  // as long as grid1 and grid2 are of the same grid type and regular level it does not matter
  //  from which the levels and indices are taken
  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  double result_ff = 0.0;
  double result_dff = 0.0;
  double result_dfdf = 0.0;
  double result_fg = 0.0;
  for (size_t i = 0; i < gridStorage.getSize(); i++) {
    sgpp::base::GridPoint& basisI = gridStorage.getPoint(i);
    size_t levelI = basisI.getLevel(0);
    size_t indexI = basisI.getIndex(0);
    for (size_t j = 0; j < gridStorage.getSize(); j++) {
      sgpp::base::GridPoint& basisJ = gridStorage.getPoint(j);
      size_t levelJ = basisJ.getLevel(0);
      size_t indexJ = basisJ.getIndex(0);
      result_ff += alpha[i] * alpha[j] *
                   ASM.univariateScalarProduct(levelI, indexI, false, levelJ, indexJ, false,
                                               pCoordinates, pWeights);
      result_dff += alpha[i] * alpha[j] *
                    ASM.univariateScalarProduct(levelI, indexI, false, levelJ, indexJ, true,
                                                pCoordinates, pWeights);
      result_dfdf += alpha[i] * alpha[j] *
                     ASM.univariateScalarProduct(levelI, indexI, true, levelJ, indexJ, true,
                                                 pCoordinates, pWeights);
      result_fg += alpha[i] * alpha2[j] *
                   ASM.univariateScalarProduct(levelI, indexI, false, levelJ, indexJ, false,
                                               pCoordinates, pWeights);
    }
  }
  double epsilon = 1e-15;
  double correctResult_ff = 1.0 / 7.0;
  double correctResult_dff = 1.0 / 2.0;
  double correctResult_dfdf = 9.0 / 5.0;
  double correctResult_fg = 0.177098574917009067047176;
  double err_ff = abs(result_ff - correctResult_ff);
  double err_dff = abs(result_dff - correctResult_dff);
  double err_dfdf = abs(result_dfdf - correctResult_dfdf);
  double err_fg = abs(result_fg - correctResult_fg);
  //  std::cout << "res ff   =" << result_ff << " err =" << err_ff << std::endl;
  //  std::cout << "res dff  =" << result_dff << " err =" << err_dff << std::endl;
  //  std::cout << "res dfdf =" << result_dfdf << " err =" << err_dfdf << std::endl;
  //  std::cout << "res fg   =" << result_fg << " err =" << err_fg << std::endl;
  BOOST_CHECK_SMALL(err_ff, epsilon);
  BOOST_CHECK_SMALL(err_dff, epsilon);
  BOOST_CHECK_SMALL(err_dfdf, epsilon);
  BOOST_CHECK_SMALL(err_fg, 1e-6);
}

BOOST_AUTO_TEST_CASE(testASMatrixEigenValuesAndVectors) {
  Eigen::MatrixXd C(3, 3);
  C(0, 0) = 1;
  C(0, 1) = 2;
  C(0, 2) = -1;
  C(1, 0) = 2;
  C(1, 1) = 5;
  C(1, 2) = 0;
  C(2, 0) = -1;
  C(2, 1) = 0;
  C(2, 2) = 3;
  auto dummyFunc = std::make_shared<sgpp::optimization::WrapperScalarFunction>(3, dummyFunction);
  sgpp::base::GridType dummyGridType = sgpp::base::GridType::NakBspline;
  size_t dummyDegree = 3;
  sgpp::optimization::ASMatrixNakBspline ASM(dummyFunc, dummyGridType, dummyDegree);
  ASM.setMatrix(C);
  ASM.evDecompositionForSymmetricMatrices();
  Eigen::VectorXd eigenvalues = ASM.getEigenvalues();
  Eigen::MatrixXd eigenvectors = ASM.getEigenvectors();
  Eigen::VectorXd ev1 = eigenvectors.col(0);
  Eigen::VectorXd ev2 = eigenvectors.col(1);
  Eigen::VectorXd ev3 = eigenvectors.col(2);

  Eigen::VectorXd Cev1 = C * ev1;
  Eigen::VectorXd e1ev1 = eigenvalues[0] * ev1;
  Eigen::VectorXd Cev2 = C * ev2;
  Eigen::VectorXd e2ev2 = eigenvalues[1] * ev2;
  Eigen::VectorXd Cev3 = C * ev3;
  Eigen::VectorXd e3ev3 = eigenvalues[2] * ev3;

  double diff1 = abs(Cev1(0) - e1ev1(0)) + abs(Cev1(1) - e1ev1(1)) + abs(Cev1(2) - e1ev1(2));
  double diff2 = abs(Cev2(0) - e2ev2(0)) + abs(Cev2(1) - e2ev2(1)) + abs(Cev2(2) - e2ev2(2));
  double diff3 = abs(Cev3(0) - e3ev3(0)) + abs(Cev3(1) - e3ev3(1)) + abs(Cev3(2) - e3ev3(2));
  double epsilon = 1e-14;
  BOOST_CHECK_SMALL(diff1, epsilon);
  BOOST_CHECK_SMALL(diff2, epsilon);
  BOOST_CHECK_SMALL(diff3, epsilon);
}

double objectiveFunctionResponseSurface(sgpp::base::DataVector v) { return v[0] + 2 * v[1]; }

BOOST_AUTO_TEST_CASE(testASResponseSurfaceNakBspline) {
  size_t numDim = 2;
  size_t degree = 3;
  size_t level = 3;
  size_t numMCPoints = 2000;

  auto objectiveFunc = std::make_shared<sgpp::optimization::WrapperScalarFunction>(
      numDim, objectiveFunctionResponseSurface);
  sgpp::base::GridType gridType = sgpp::base::GridType::NakBsplineBoundary;
  sgpp::optimization::ASMatrixNakBspline ASM(objectiveFunc, gridType, degree);
  ASM.buildRegularInterpolant(level);
  ASM.createMatrixMonteCarlo(numMCPoints);
  ASM.evDecompositionForSymmetricMatrices();
  Eigen::VectorXd eigenvalues = ASM.getEigenvalues();
  Eigen::MatrixXd eigenvectors = ASM.getEigenvectors();
  // active subspace specifier
  size_t n = 1;
  Eigen::MatrixXd W1 = ASM.getTransformationMatrix(n);

  sgpp::optimization::ASResponseSurfaceNakBspline responseSurf(numDim, W1, gridType, degree);
  sgpp::base::DataMatrix evaluationPoints = ASM.getEvaluationPoints();
  sgpp::base::DataVector functionValues = ASM.getFunctionValues();
  size_t responseLevel = 3;
  responseSurf.createRegularReducedSurfaceFromDetectionPoints(evaluationPoints, functionValues,
                                                              responseLevel);
  sgpp::base::DataVector v(numDim, 0.3371);
  double objectiveFunctionEval = objectiveFunctionResponseSurface(v);
  double responseSurfEval = responseSurf.eval(v);
  double epsilon = 1e-15;
  BOOST_CHECK_SMALL(fabs(objectiveFunctionEval - responseSurfEval), epsilon);
}

/* ToDo (rehmemk)
 * -Test response surface
 */

BOOST_AUTO_TEST_SUITE_END()
