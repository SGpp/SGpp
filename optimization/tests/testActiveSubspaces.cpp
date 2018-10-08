// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp/optimization/activeSubspaces/GaussQuadrature.hpp>
#include <sgpp/optimization/tools/Printer.hpp>

#include <functional>

double oneFunction(double x) { return 1; }

double xFunction(double x) { return x; }

double x6Function(double x) { return x * x * x * x * x * x; }

BOOST_AUTO_TEST_SUITE(testActiveSubspaces)

BOOST_AUTO_TEST_CASE(TestQuad) {
  double epsilon = 1e-15;

  std::function<double(double)> func = oneFunction;
  size_t quadOrder = 1;
  double correctRes = 5;
  double quadResult = sgpp::optimization::quad(func, -3, 2, quadOrder);
  double err = fabs(correctRes - quadResult);
  std::cout << err << std::endl;
  BOOST_CHECK_SMALL(err, epsilon);

  func = xFunction;
  quadOrder = 1;
  correctRes = 50;
  quadResult = sgpp::optimization::quad(func, 0, 10, quadOrder);
  err = fabs(correctRes - quadResult);
  std::cout << err << std::endl;
  BOOST_CHECK_SMALL(err, epsilon);

  func = x6Function;
  quadOrder = 4;
  correctRes = 1.0 / 7.0;
  quadResult = sgpp::optimization::quad(func, 0, 1, quadOrder);
  err = fabs(correctRes - quadResult);
  std::cout << err << std::endl;
  BOOST_CHECK_SMALL(err, epsilon);
}

BOOST_AUTO_TEST_CASE(TestThat) { sgpp::optimization::Printer::getInstance().setVerbosity(-1); }

BOOST_AUTO_TEST_SUITE_END()
