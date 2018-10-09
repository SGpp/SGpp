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

double objectiveFunction(double x) { return x * x * x; }
double dummyFunction(sgpp::base::DataVector v) { return 777; }

BOOST_AUTO_TEST_CASE(testASMatrixNakBsplineBoundaryScalarProduct) {
  // interpolate function f, calculate \int f^2 dx with the scalar product routines of
  // ASMatrixNakBspline
  sgpp::optimization::Printer::getInstance().setVerbosity(-1);

  size_t numDim = 1;
  size_t degree = 3;
  size_t level = 3;

  auto grid = std::make_shared<sgpp::base::NakBsplineBoundaryGrid>(numDim, degree);
  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  grid->getGenerator().regular(level);
  sgpp::base::DataVector functionValues(gridStorage.getSize());
  for (size_t i = 0; i < gridStorage.getSize(); i++) {
    sgpp::base::GridPoint& gp = gridStorage.getPoint(i);
    functionValues[i] = objectiveFunction(gp.getStandardCoordinate(0));
  }

  // solve linear system
  sgpp::base::DataVector alpha(functionValues.getSize());
  sgpp::optimization::HierarchisationSLE hierSLE(*grid);
  sgpp::optimization::sle_solver::Armadillo sleSolver;
  if (!sleSolver.solve(hierSLE, functionValues, alpha)) {
    std::cout << "TestASMatrixNakBsplineScalarProduct: Solving failed.\n";
    return;
  }

  sgpp::optimization::WrapperScalarFunction dummyFunc(numDim, dummyFunction);
  sgpp::base::GridType gridType = sgpp::base::GridType::NakBsplineBoundary;
  sgpp::optimization::ASMatrixNakBspline ASM(dummyFunc, gridType, degree);

  double result = 0.0;
  for (size_t i = 0; i < gridStorage.getSize(); i++) {
    sgpp::base::GridPoint& basisI = gridStorage.getPoint(i);
    size_t levelI = basisI.getLevel(0);
    size_t indexI = basisI.getIndex(0);
    for (size_t j = 0; j < gridStorage.getSize(); j++) {
      sgpp::base::GridPoint& basisJ = gridStorage.getPoint(j);
      size_t levelJ = basisJ.getLevel(0);
      size_t indexJ = basisJ.getIndex(0);
      result += alpha[i] * alpha[j] *
                ASM.univariateScalarProduct(levelI, indexI, false, levelJ, indexJ, false);
    }
  }
  double epsilon = 1e-15;
  double correctResult = 1.0 / 7.0;
  double err = abs(result - correctResult);
  std::cout << "res=" << result << " err=" << err << std::endl;
  BOOST_CHECK_SMALL(err, epsilon);
}

BOOST_AUTO_TEST_SUITE_END()
