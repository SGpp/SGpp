// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <sgpp/base/function/scalar/InterpolantScalarFunction.hpp>
#include <sgpp/base/function/scalar/WrapperScalarFunction.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/type/NakBsplineExtendedGrid.hpp>
#include <sgpp/base/tools/Printer.hpp>
#include <sgpp/base/tools/RandomNumberGenerator.hpp>
//#include <sgpp/base/tools/sle/solver/Armadillo.hpp>
#include <sgpp/base/tools/sle/solver/Auto.hpp>
#include <sgpp/base/tools/sle/system/HierarchisationSLE.hpp>

BOOST_AUTO_TEST_SUITE(TestBasis)

// test functions
double func1(sgpp::base::DataVector v) { return 1; }
double funcx(sgpp::base::DataVector v) { return v[0]; }
double funcx2(sgpp::base::DataVector v) { return v[0] * v[0]; }
double funcx3(sgpp::base::DataVector v) { return v[0] * v[0] * v[0]; }
double funcx4(sgpp::base::DataVector v) { return v[0] * v[0] * v[0] * v[0]; }
double funcx5(sgpp::base::DataVector v) { return v[0] * v[0] * v[0] * v[0] * v[0]; }

// auxiliary function. Interpolates the objective function func on a regular grid of level level
// with a grid of GridType and degree and return the l2 error caluclated on numMCPoints Monte Carlo
// points.
// returns grid and alpha
void interpolate(size_t degree, size_t level, sgpp::base::WrapperScalarFunction func,
                 sgpp::base::GridType gridType, std::shared_ptr<sgpp::base::Grid>& grid,
                 sgpp::base::DataVector& alpha) {
  size_t dim = 1;
  if (gridType == sgpp::base::GridType::NakBsplineExtended) {
    grid.reset(new sgpp::base::NakBsplineExtendedGrid(dim, degree));
  } else {
    std::cerr << "test_Basis: interpolate: gridType not supported";
  }
  grid->getGenerator().regular(level);
  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  sgpp::base::DataVector f_values(gridStorage.getSize());
  for (size_t i = 0; i < gridStorage.getSize(); i++) {
    sgpp::base::GridPoint& gp = gridStorage.getPoint(i);
    sgpp::base::DataVector gridPointVector(gridStorage.getDimension());
    gp.getStandardCoordinates(gridPointVector);
    f_values[i] = func.eval(gridPointVector);
  }
  alpha.resizeZero(gridStorage.getSize());
  sgpp::base::HierarchisationSLE hierSLE(*grid);
  // sgpp::base::sle_solver::Armadillo sleSolver;
  sgpp::base::sle_solver::Auto sleSolver;
  if (!sleSolver.solve(hierSLE, f_values, alpha)) {
    std::cout << "Solving failed.\n";
    return;
  }
}
double l2Error(std::shared_ptr<sgpp::base::Grid> grid, sgpp::base::DataVector alpha,
               sgpp::base::WrapperScalarFunction func, size_t numMCPoints) {
  sgpp::base::InterpolantScalarFunction I(*grid, alpha);
  double l2Err = 0.0;
  sgpp::base::DataVector randomVector(func.getNumberOfParameters());
  for (size_t i = 0; i < numMCPoints; i++) {
    sgpp::base::RandomNumberGenerator::getInstance().getUniformRV(randomVector, 0.0, 1.0);
    double evalInterpolant = I.eval(randomVector);
    double evalObjectiveFunc = func.eval(randomVector);
    l2Err += std::pow(evalInterpolant - evalObjectiveFunc, 2.0);
  }
  l2Err = sqrt(l2Err / static_cast<double>(numMCPoints));
  return l2Err;
}

double integrate(std::shared_ptr<sgpp::base::Grid> grid, sgpp::base::DataVector alpha,
                 size_t degree) {
  std::shared_ptr<sgpp::base::SBasis> basis;
  if (grid->getType() == sgpp::base::GridType::NakBsplineExtended) {
    basis.reset(new sgpp::base::SNakBsplineExtendedBase(degree));
  } else {
    std::cerr << "test_Basis: integrate: gridType not supported";
  }
  double integral = 0;
  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  for (size_t j = 0; j < gridStorage.getSize(); j++) {
    sgpp::base::GridPoint& gpBasis = gridStorage.getPoint(j);
    double basisInt = 1;
    for (size_t t = 0; t < gridStorage.getDimension(); t++) {
      double basisInt1D = basis->getIntegral(gpBasis.getLevel(t), gpBasis.getIndex(t));
      basisInt *= basisInt1D;
    }
    integral += alpha[j] * basisInt;
  }
  return integral;
}

void testInterpolationAndIntegration(sgpp::base::GridType gridType) {
  // test if Basis of degree p in {1,3,5} represents polynomials of degree p exact
  sgpp::base::Printer::getInstance().setVerbosity(-1);
  size_t numMCPoints = 1000;
  std::shared_ptr<sgpp::base::Grid> grid;
  sgpp::base::DataVector alpha;
  sgpp::base::WrapperScalarFunction Func1(1, func1);
  sgpp::base::WrapperScalarFunction Funcx(1, funcx);
  sgpp::base::WrapperScalarFunction Funcx2(1, funcx2);
  sgpp::base::WrapperScalarFunction Funcx3(1, funcx3);
  sgpp::base::WrapperScalarFunction Funcx4(1, funcx4);
  sgpp::base::WrapperScalarFunction Funcx5(1, funcx5);

  interpolate(1, 1, Func1, gridType, grid, alpha);
  double l2error11 = l2Error(grid, alpha, Func1, numMCPoints);
  double interror11 = abs(integrate(grid, alpha, 1) - 1);
  interpolate(1, 2, Funcx, gridType, grid, alpha);
  double l2error1x = l2Error(grid, alpha, Funcx, numMCPoints);
  double interror1x = abs(integrate(grid, alpha, 1) - 0.5);

  interpolate(3, 1, Func1, gridType, grid, alpha);
  double l2error31 = l2Error(grid, alpha, Func1, numMCPoints);
  double interror31 = abs(integrate(grid, alpha, 3) - 1);
  interpolate(3, 2, Funcx, gridType, grid, alpha);
  double l2error3x = l2Error(grid, alpha, Funcx, numMCPoints);
  double interror3x = abs(integrate(grid, alpha, 3) - 0.5);
  interpolate(3, 2, Funcx2, gridType, grid, alpha);
  double l2error3x2 = l2Error(grid, alpha, Funcx2, numMCPoints);
  double interror3x2 = abs(integrate(grid, alpha, 3) - 1.0 / 3.0);
  interpolate(3, 3, Funcx3, gridType, grid, alpha);
  double l2error3x3 = l2Error(grid, alpha, Funcx3, numMCPoints);
  double interror3x3 = abs(integrate(grid, alpha, 3) - 0.25);

  interpolate(5, 1, Func1, gridType, grid, alpha);
  double l2error51 = l2Error(grid, alpha, Func1, numMCPoints);
  double interror51 = abs(integrate(grid, alpha, 5) - 1);
  interpolate(5, 2, Funcx, gridType, grid, alpha);
  double l2error5x = l2Error(grid, alpha, Funcx, numMCPoints);
  double interror5x = abs(integrate(grid, alpha, 5) - 0.5);
  interpolate(5, 2, Funcx2, gridType, grid, alpha);
  double l2error5x2 = l2Error(grid, alpha, Funcx2, numMCPoints);
  double interror5x2 = abs(integrate(grid, alpha, 5) - 1.0 / 3.0);
  interpolate(5, 3, Funcx3, gridType, grid, alpha);
  double l2error5x3 = l2Error(grid, alpha, Funcx3, numMCPoints);
  double interror5x3 = abs(integrate(grid, alpha, 5) - 0.25);
  interpolate(5, 3, Funcx4, gridType, grid, alpha);
  double l2error5x4 = l2Error(grid, alpha, Funcx4, numMCPoints);
  double interror5x4 = abs(integrate(grid, alpha, 5) - 1.0 / 5.0);
  interpolate(5, 3, Funcx5, gridType, grid, alpha);
  double l2error5x5 = l2Error(grid, alpha, Funcx5, numMCPoints);
  double interror5x5 = abs(integrate(grid, alpha, 5) - 1.0 / 6.0);

  //  std::cout << l2error11 << " " << l2error1x << "\n";
  //  std::cout << l2error31 << " " << l2error3x << " " << l2error3x2 << " " << l2error3x3 << "\n";
  //  std::cout << l2error51 << " " << l2error5x << " " << l2error5x2 << " " << l2error5x3 << " "
  //            << l2error5x4 << " " << l2error5x5 << "\n";

  //  std::cout << interror11 << " " << interror1x << "\n";
  //  std::cout << interror31 << " " << interror3x << " " << interror3x2 << " " << interror3x3 <<
  //  "\n"; std::cout << interror51 << " " << interror5x << " " << interror5x2 << " " << interror5x3
  //  << " "
  //            << interror5x4 << " " << interror5x5 << "\n";

  double tolerance = 1e-15;
  BOOST_CHECK_SMALL(l2error11, tolerance);
  BOOST_CHECK_SMALL(l2error1x, tolerance);

  BOOST_CHECK_SMALL(l2error31, tolerance);
  BOOST_CHECK_SMALL(l2error3x, tolerance);
  BOOST_CHECK_SMALL(l2error3x2, tolerance);
  BOOST_CHECK_SMALL(l2error3x3, tolerance);

  BOOST_CHECK_SMALL(l2error51, tolerance);
  BOOST_CHECK_SMALL(l2error5x, tolerance);
  BOOST_CHECK_SMALL(l2error5x2, tolerance);
  BOOST_CHECK_SMALL(l2error5x3, tolerance);
  BOOST_CHECK_SMALL(l2error5x4, tolerance);
  BOOST_CHECK_SMALL(l2error5x5, tolerance);

  BOOST_CHECK_SMALL(interror11, tolerance);
  BOOST_CHECK_SMALL(interror1x, tolerance);

  BOOST_CHECK_SMALL(interror31, tolerance);
  BOOST_CHECK_SMALL(interror3x, tolerance);
  BOOST_CHECK_SMALL(interror3x2, tolerance);
  BOOST_CHECK_SMALL(interror3x3, tolerance);

  BOOST_CHECK_SMALL(interror51, tolerance);
  BOOST_CHECK_SMALL(interror5x, tolerance);
  BOOST_CHECK_SMALL(interror5x2, tolerance);
  BOOST_CHECK_SMALL(interror5x3, tolerance);
  BOOST_CHECK_SMALL(interror5x4, tolerance);
  BOOST_CHECK_SMALL(interror5x5, tolerance);
}

BOOST_AUTO_TEST_CASE(testNakBsplineExtendedBasis) {
  sgpp::base::GridType gridType = sgpp::base::GridType::NakBsplineExtended;
  testInterpolationAndIntegration(gridType);
}

BOOST_AUTO_TEST_SUITE_END()
