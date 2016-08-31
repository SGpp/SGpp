// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp/combigrid/grid/ordering/ExponentialChebyshevPointOrdering.hpp>
#include <sgpp/combigrid/grid/ordering/ExponentialLevelorderPointOrdering.hpp>
#include <sgpp/combigrid/grid/ordering/ExponentialNoBoundaryPointOrdering.hpp>
#include <vector>

BOOST_AUTO_TEST_CASE(testExpLevelorderPointOrdering) {
  sgpp::combigrid::ExponentialLevelorderPointOrdering ordering;

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

  std::vector<double> points;

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

BOOST_AUTO_TEST_CASE(testExpChebyshevPointOrdering) {
  sgpp::combigrid::ExponentialChebyshevPointOrdering ordering;

  BOOST_CHECK_EQUAL(ordering.numPoints(0), 1);
  BOOST_CHECK_EQUAL(ordering.numPoints(1), 3);
  BOOST_CHECK_EQUAL(ordering.numPoints(2), 9);
  BOOST_CHECK_EQUAL(ordering.numPoints(3), 27);

  BOOST_CHECK_EQUAL(ordering.convertIndex(2, 9, 0), 4);
  BOOST_CHECK_EQUAL(ordering.convertIndex(2, 9, 1), 1);
  BOOST_CHECK_EQUAL(ordering.convertIndex(2, 9, 2), 7);
  BOOST_CHECK_EQUAL(ordering.convertIndex(2, 9, 3), 0);
  BOOST_CHECK_EQUAL(ordering.convertIndex(2, 9, 4), 2);
  BOOST_CHECK_EQUAL(ordering.convertIndex(2, 9, 5), 3);
  BOOST_CHECK_EQUAL(ordering.convertIndex(2, 9, 6), 5);
  BOOST_CHECK_EQUAL(ordering.convertIndex(2, 9, 7), 6);
  BOOST_CHECK_EQUAL(ordering.convertIndex(2, 9, 8), 8);

  std::vector<double> points;

  auto it = ordering.getSortedPermutationIterator(2, points, 9);

  BOOST_CHECK_EQUAL(it->value(), 3);
  it->moveToNext();
  BOOST_CHECK_EQUAL(it->value(), 1);
  it->moveToNext();
  BOOST_CHECK_EQUAL(it->value(), 4);
  it->moveToNext();
  BOOST_CHECK_EQUAL(it->value(), 5);
  it->moveToNext();
  BOOST_CHECK_EQUAL(it->value(), 0);
  it->moveToNext();
  BOOST_CHECK_EQUAL(it->value(), 6);
  it->moveToNext();
  BOOST_CHECK_EQUAL(it->value(), 7);
  it->moveToNext();
  BOOST_CHECK_EQUAL(it->value(), 2);
  it->moveToNext();
  BOOST_CHECK_EQUAL(it->value(), 8);
  it->reset();
  BOOST_CHECK_EQUAL(it->value(), 3);
}

BOOST_AUTO_TEST_CASE(testExpNoBoundaryPointOrdering) {
  sgpp::combigrid::ExponentialNoBoundaryPointOrdering ordering;

  BOOST_CHECK_EQUAL(ordering.numPoints(0), 1);
  BOOST_CHECK_EQUAL(ordering.numPoints(1), 3);
  BOOST_CHECK_EQUAL(ordering.numPoints(2), 7);
  BOOST_CHECK_EQUAL(ordering.numPoints(3), 15);

  BOOST_CHECK_EQUAL(ordering.convertIndex(2, 7, 0), 3);
  BOOST_CHECK_EQUAL(ordering.convertIndex(2, 7, 1), 1);
  BOOST_CHECK_EQUAL(ordering.convertIndex(2, 7, 2), 5);
  BOOST_CHECK_EQUAL(ordering.convertIndex(2, 7, 3), 0);
  BOOST_CHECK_EQUAL(ordering.convertIndex(2, 7, 4), 2);
  BOOST_CHECK_EQUAL(ordering.convertIndex(2, 7, 5), 4);
  BOOST_CHECK_EQUAL(ordering.convertIndex(2, 7, 6), 6);

  std::vector<double> points;

  auto it = ordering.getSortedPermutationIterator(2, points, 7);

  BOOST_CHECK_EQUAL(it->value(), 3);
  it->moveToNext();
  BOOST_CHECK_EQUAL(it->value(), 1);
  it->moveToNext();
  BOOST_CHECK_EQUAL(it->value(), 4);
  it->moveToNext();
  BOOST_CHECK_EQUAL(it->value(), 0);
  it->moveToNext();
  BOOST_CHECK_EQUAL(it->value(), 5);
  it->moveToNext();
  BOOST_CHECK_EQUAL(it->value(), 2);
  it->moveToNext();
  BOOST_CHECK_EQUAL(it->value(), 6);
  it->reset();
  BOOST_CHECK_EQUAL(it->value(), 3);
}
