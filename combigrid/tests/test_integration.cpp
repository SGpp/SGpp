// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp/combigrid/operation/onedim/QuadratureEvaluator.hpp>
#include <sgpp/combigrid/definitions.hpp>
#include <sgpp/globaldef.hpp>

#include <vector>

const double tolerance = 1e-12;

BOOST_AUTO_TEST_CASE(testIntegration) {
  sgpp::combigrid::QuadratureEvaluator eval;

  std::vector<double> points;
  std::vector<double> weights;

  // 1 Point

  points.push_back(0.0);

  eval.setGridPoints(points);

  BOOST_CHECK_EQUAL(eval.getBasisCoefficients().size(), 1);
  BOOST_CHECK_CLOSE(eval.getBasisCoefficients()[0].getValue(), 1.0, tolerance);

  // 2 Points

  points.push_back(1.0);

  eval.setGridPoints(points);

  BOOST_CHECK_EQUAL(eval.getBasisCoefficients().size(), 2);

  BOOST_CHECK_CLOSE(eval.getBasisCoefficients()[0].getValue(), 0.5, tolerance);
  BOOST_CHECK_CLOSE(eval.getBasisCoefficients()[1].getValue(), 0.5, tolerance);
}

void testIntegrationWithPolynom() {
  sgpp::combigrid::QuadratureEvaluator eval;

  std::vector<double> points;
  std::vector<double> weights;

  points.push_back(0.0);
  points.push_back(0.5);
  points.push_back(1.0);

  eval.setGridPoints(points);

  // p(x) = x^2

  std::vector<double> function_points;
  function_points.push_back(0.0);
  function_points.push_back(0.25);
  function_points.push_back(1.0);

  // check coefficients
  BOOST_CHECK_EQUAL(eval.getBasisCoefficients().size(), 3);
  BOOST_CHECK_CLOSE(eval.getBasisCoefficients()[0].getValue(), 1.0 / 6, tolerance);
  BOOST_CHECK_CLOSE(eval.getBasisCoefficients()[1].getValue(), 2.0 / 3, tolerance);
  BOOST_CHECK_CLOSE(eval.getBasisCoefficients()[2].getValue(), 1.0 / 6, tolerance);

  double accurate_solution = 1.0 / 3;

  BOOST_CHECK_CLOSE(eval.eval(function_points).getValue(), accurate_solution, tolerance);

  // try a more complex polynomial
  // p(x) = 6 * x^5 - 5 * x^4 + 4 * x^3 - 3 * x^2 + 2 * x - 42
  // the accurate integral is x * (x^5 - x^4 + x^3 - x^2 + x - 42)
  // so in [0,1], it is -41
  accurate_solution = -41.0;

  // we need to give d + 1 points, here d is 5
  points.clear();
  function_points.clear();
  points.push_back(0.0);
  function_points.push_back(-42.0);
  points.push_back(0.2);
  function_points.push_back(-41.69408);
  points.push_back(0.4);
  function_points.push_back(-41.49056);
  points.push_back(0.6);
  function_points.push_back(-41.19744);
  points.push_back(0.8);
  function_points.push_back(-40.35392);
  points.push_back(1.0);
  function_points.push_back(-38.0);

  eval.setGridPoints(points);

  BOOST_CHECK_CLOSE(eval.eval(function_points).getValue(), accurate_solution, 0.2);
}

BOOST_AUTO_TEST_CASE(testQuadrature) { testIntegrationWithPolynom(); }
