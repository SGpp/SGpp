/*
 * test_leja_points.cpp
 *
 *  Created on: Jan 25, 2016
 *      Author: liedtkjn
 */

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp/combigrid/grid/distribution/LejaPointDistribution.hpp>
#include <cmath>

// in %
SGPP::float_t tolerance = 0.01;

void testLejaPoints() {
	std::function<SGPP::float_t(SGPP::float_t)> weight_func =
			[](SGPP::float_t d) {return 1;};
	SGPP::float_t start = 0.5;
	SGPP::float_t lower_bound = 0;
	SGPP::float_t upper_bound = 1;
	int number = 8;

	std::vector<SGPP::float_t> leja_points;

	// start point
	leja_points.push_back(start);
	auto points = leja_points;

	// calc leja points
	SGPP::combigrid::calc_leja_points(leja_points, points, number, lower_bound,
			upper_bound, weight_func);

	// correct solution:
	std::vector<SGPP::float_t> correct_leja_points;
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
		BOOST_CHECK_CLOSE(leja_points.at(i), correct_leja_points.at(i),
				tolerance);
	}
}

void testLejaDistributionStartingPoint() {
	SGPP::combigrid::LejaPointDistribution leja;

	SGPP::float_t first_point = leja.compute(0, 0);

	BOOST_CHECK_CLOSE(first_point, 0.5, tolerance);
}

/*
 * the starting point for this function is 1.0
 */
SGPP::float_t quadratic_func(SGPP::float_t x) {
	return 42 * x * x;
}

void testOtherStartingPoint() {
	SGPP::combigrid::LejaPointDistribution leja(quadratic_func);

	// check if first point is correct, it should be 1.0
	BOOST_CHECK_CLOSE(leja.compute(0, 0), 1.0, tolerance);
}

SGPP::float_t sinusweight(SGPP::float_t x) {
	return std::sin(x * 3.14159265358979323846);
}

void testLejaSinusWithNormalDistribution() {
	SGPP::combigrid::LejaPointDistribution leja(sinusweight);

	// check if first point is correct, it should be 0.5
	BOOST_CHECK_CLOSE(leja.compute(0, 0), 0.5, tolerance);

	// check the next few points
	BOOST_CHECK_CLOSE(leja.compute(1, 1), 0.22614641471525609, tolerance);
	BOOST_CHECK_CLOSE(leja.compute(1, 2), 0.81695850918945612, tolerance);
	BOOST_CHECK_CLOSE(leja.compute(1, 3), 0.088181946121442964, tolerance);
	BOOST_CHECK_CLOSE(leja.compute(1, 4), 0.92906611726964328, tolerance);
}

/**
 * Warning: This tests are made for a epsilon of 0.00001 in the leja class (for the optimizer)!
 */
BOOST_AUTO_TEST_CASE(testLeja) {
	testLejaPoints();
	testLejaDistributionStartingPoint();
	// check what happens if the starting point isn't 0.5
	testOtherStartingPoint();
	// more complex test case
	testLejaSinusWithNormalDistribution();
}
