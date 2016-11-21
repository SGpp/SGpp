// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp/combigrid/SingleFunction.hpp>
#include <sgpp/combigrid/grid/distribution/LejaPointDistribution.hpp>

#include <cmath>
#include <vector>

// in %
double tolerance = 0.01;

void testLejaPoints() {
  std::function<double(double)> weight_func = [](double d) { return 1; };
  double start = 0.5;
  double lower_bound = 0;
  double upper_bound = 1;
  int number = 8;

  std::vector<double> leja_points;

  // start point
  leja_points.push_back(start);
  auto points = leja_points;

  // calc leja points
  sgpp::combigrid::LejaPointDistribution::calc_leja_points(leja_points, points, number, lower_bound,
                                                           upper_bound, weight_func);

  // correct solution:
  std::vector<double> correct_leja_points;
  correct_leja_points.push_back(1.0428174199820221e-05);
  correct_leja_points.push_back(0.065002567709224635);
  correct_leja_points.push_back(0.17065421163314368);
  correct_leja_points.push_back(0.34719663104769344);
  correct_leja_points.push_back(0.5);
  correct_leja_points.push_back(0.66085325446028498);
  correct_leja_points.push_back(0.78867243975866386);
  correct_leja_points.push_back(0.91962400339636541);
  correct_leja_points.push_back(0.99998958672635252);

  // test number of points
  BOOST_CHECK_EQUAL(leja_points.size(), correct_leja_points.size());

  for (size_t i = 0; i < correct_leja_points.size(); ++i) {
    BOOST_CHECK_CLOSE(leja_points.at(i), correct_leja_points.at(i), tolerance);
  }
}

void testLejaDistributionStartingPoint() {
  sgpp::combigrid::LejaPointDistribution leja;

  double first_point = leja.compute(0, 0);

  BOOST_CHECK_CLOSE(first_point, 0.5, tolerance);
}

/*
 * the starting point for this function is 1.0
 */
double quadratic_func(double x) { return 42 * x * x; }

void testOtherStartingPoint() {
  auto w = sgpp::combigrid::SingleFunction(quadratic_func);
  sgpp::combigrid::LejaPointDistribution leja(w);

  // check if first point is correct, it should be 1.0
  BOOST_CHECK_CLOSE(leja.compute(0, 0), 1.0, tolerance);
}

double sinusweight(double x) { return std::sin(x * 3.14159265358979323846); }

void testLejaSinusWithNormalDistribution() {
  auto w = sgpp::combigrid::SingleFunction(sinusweight);
  sgpp::combigrid::LejaPointDistribution leja(w);

  // check if first point is correct, it should be 0.5
  BOOST_CHECK_CLOSE(leja.compute(0, 0), 0.5, tolerance);

  // check the next few points
  BOOST_CHECK_CLOSE(leja.compute(1, 1), 1 - 0.22614641471525609, tolerance);
  BOOST_CHECK_CLOSE(leja.compute(1, 2), 1 - 0.81695850918945612, tolerance);
  BOOST_CHECK_CLOSE(leja.compute(1, 3), 1 - 0.088181946121442964, tolerance);
  BOOST_CHECK_CLOSE(leja.compute(1, 4), 1 - 0.92906611726964328, tolerance);
}

/**
 * Warning: This tests are made for a epsilon of 0.00001 in the leja class (for the optimizer)!
 */
BOOST_AUTO_TEST_CASE(testLeja) {
  //  testLejaPoints();
  testLejaDistributionStartingPoint();
  // check what happens if the starting point isn't 0.5
  testOtherStartingPoint();
  // more complex test case
  testLejaSinusWithNormalDistribution();
}
