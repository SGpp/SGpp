/*
 * test_point_ordering.cpp
 *
 *  Created on: 19.12.2015
 *      Author: david
 */

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <sgpp/combigrid/grid/ordering/ExponentialLevelorderPointOrdering.hpp>

using namespace SGPP::combigrid;

BOOST_AUTO_TEST_CASE(testExpLevelorderPointOrdering) {
	ExponentialLevelorderPointOrdering ordering;

	BOOST_CHECK_EQUAL(ordering.numPoints(0), 1);
	BOOST_CHECK_EQUAL(ordering.numPoints(1), 3);
	BOOST_CHECK_EQUAL(ordering.numPoints(2), 5);
	BOOST_CHECK_EQUAL(ordering.numPoints(3), 9);

	BOOST_CHECK_EQUAL(ordering.convertIndex(3, 9, 0), 4);
	BOOST_CHECK_EQUAL(ordering.convertIndex(3, 9, 1), 0);
	BOOST_CHECK_EQUAL(ordering.convertIndex(3, 9, 2), 8);
	BOOST_CHECK_EQUAL(ordering.convertIndex(3, 9, 3), 2);
	BOOST_CHECK_EQUAL(ordering.convertIndex(3, 9, 4), 6);
	BOOST_CHECK_EQUAL(ordering.convertIndex(3, 9, 5), 1);
	BOOST_CHECK_EQUAL(ordering.convertIndex(3, 9, 6), 3);
	BOOST_CHECK_EQUAL(ordering.convertIndex(3, 9, 7), 5);
	BOOST_CHECK_EQUAL(ordering.convertIndex(3, 9, 8), 7);

	std::vector<SGPP::float_t> points;

	auto it = ordering.getSortedPermutationIterator(3, points, 9);

	BOOST_CHECK_EQUAL(it->value(), 1);
	it->moveToNext();
	BOOST_CHECK_EQUAL(it->value(), 5);
	it->moveToNext();
	BOOST_CHECK_EQUAL(it->value(), 3);
	it->moveToNext();
	BOOST_CHECK_EQUAL(it->value(), 6);
	it->moveToNext();
	BOOST_CHECK_EQUAL(it->value(), 0);
	it->moveToNext();
	BOOST_CHECK_EQUAL(it->value(), 7);
	it->moveToNext();
	BOOST_CHECK_EQUAL(it->value(), 4);
	it->moveToNext();
	BOOST_CHECK_EQUAL(it->value(), 8);
	it->moveToNext();
	BOOST_CHECK_EQUAL(it->value(), 2);
	it->reset();
	BOOST_CHECK_EQUAL(it->value(), 1);
}


