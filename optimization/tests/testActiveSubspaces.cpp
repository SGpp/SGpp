// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp/base/grid/type/NakBsplineBoundaryGrid.hpp>
#include <sgpp/base/grid/type/NakBsplineModifiedGrid.hpp>
#include <sgpp/optimization/activeSubspaces/GaussQuadrature.hpp>
#include <sgpp/optimization/function/scalar/InterpolantScalarFunction.hpp>
#include <sgpp/optimization/function/scalar/WrapperScalarFunction.hpp>
#include <sgpp/optimization/tools/Printer.hpp>
#include "../src/sgpp/optimization/activeSubspaces/ASMatrixNakBspline.hpp"

#include <functional>

double oneFunction(double x) { return 1; }

double xFunction(double x) { return x; }

double x6Function(double x) { return x * x * x * x * x * x; }

BOOST_AUTO_TEST_SUITE(testActiveSubspaces)

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
    sgpp::optimization::WrapperScalarFunction objectiveFunction, sgpp::base::GridType gridType,
    size_t degree, size_t level, sgpp::base::Grid* grid) {
  sgpp::optimization::Printer::getInstance().setVerbosity(-1);
  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  grid->getGenerator().regular(level);
  sgpp::base::DataVector functionValues(gridStorage.getSize());
  for (size_t i = 0; i < gridStorage.getSize(); i++) {
    sgpp::base::GridPoint& gp = gridStorage.getPoint(i);
    functionValues[i] =
        objectiveFunction.eval(sgpp::base::DataVector(1, gp.getStandardCoordinate(0)));
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

double objectiveFunction(sgpp::base::DataVector v) { return v[0] * v[0] * v[0]; }
double dummyFunction(sgpp::base::DataVector v) { return 777; }

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
  if (gridType == sgpp::base::GridType::NakBspline) {
    grid = new sgpp::base::NakBsplineGrid(numDim, degree);
  } else if (gridType == sgpp::base::GridType::NakBsplineModified) {
    grid = new sgpp::base::NakBsplineModifiedGrid(numDim, degree);
  } else if (gridType == sgpp::base::GridType::NakBsplineBoundary) {
    grid = new sgpp::base::NakBsplineBoundaryGrid(numDim, degree);
  } else {
    throw sgpp::base::generation_exception(
        "testASMatrixNakBsplineBoundaryScalarProducte: gridType not supported.");
  }
  sgpp::optimization::WrapperScalarFunction objectiveFunc(numDim, objectiveFunction);
  sgpp::base::DataVector alpha = interpolateRegular1D(objectiveFunc, gridType, degree, level, grid);
  sgpp::optimization::ASMatrixNakBspline ASM(objectiveFunc, gridType, degree);

  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  double result_ff = 0.0;
  double result_dff = 0.0;
  double result_dfdf = 0.0;
  for (size_t i = 0; i < gridStorage.getSize(); i++) {
    sgpp::base::GridPoint& basisI = gridStorage.getPoint(i);
    size_t levelI = basisI.getLevel(0);
    size_t indexI = basisI.getIndex(0);
    for (size_t j = 0; j < gridStorage.getSize(); j++) {
      sgpp::base::GridPoint& basisJ = gridStorage.getPoint(j);
      size_t levelJ = basisJ.getLevel(0);
      size_t indexJ = basisJ.getIndex(0);
      result_ff += alpha[i] * alpha[j] *
                   ASM.univariateScalarProduct(levelI, indexI, false, levelJ, indexJ, false);
      result_dff += alpha[i] * alpha[j] *
                    ASM.univariateScalarProduct(levelI, indexI, false, levelJ, indexJ, true);
      result_dfdf += alpha[i] * alpha[j] *
                     ASM.univariateScalarProduct(levelI, indexI, true, levelJ, indexJ, true);
    }
  }
  double epsilon = 1e-15;
  double correctResult_ff = 1.0 / 7.0;
  double correctResult_dff = 1.0 / 2.0;
  double correctResult_dfdf = 9.0 / 5.0;
  double err_ff = abs(result_ff - correctResult_ff);
  double err_dff = abs(result_dff - correctResult_dff);
  double err_dfdf = abs(result_dfdf - correctResult_dfdf);
  //  std::cout << "res ff   =" << result_ff << " err =" << err_ff << std::endl;
  //  std::cout << "res dff  =" << result_dff << " err =" << err_dff << std::endl;
  //  std::cout << "res dfdf =" << result_dfdf << " err =" << err_dfdf << std::endl;
  BOOST_CHECK_SMALL(err_ff, epsilon);
  BOOST_CHECK_SMALL(err_dff, epsilon);
  BOOST_CHECK_SMALL(err_dfdf, epsilon);
}

BOOST_AUTO_TEST_SUITE_END()
