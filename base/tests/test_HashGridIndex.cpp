// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <sgpp/base/grid/storage/hashmap/HashGridIndex.hpp>

using sgpp::base::HashGridIndex;

BOOST_AUTO_TEST_SUITE(TestSHashGridIndex)

BOOST_AUTO_TEST_CASE(testConstructor) {
  HashGridIndex s(2);
  s.set(0, 1, 1);
  s.set(0, 2, 3);

  HashGridIndex s2(s);

  BOOST_CHECK_EQUAL(s.getLevel(0), s2.getLevel(0));
  BOOST_CHECK_EQUAL(s.getIndex(0), s2.getIndex(0));
}

BOOST_AUTO_TEST_CASE(testAssign) {
  HashGridIndex s(2);
  s.set(0, 2, 1);
  s.set(1, 2, 3);

  HashGridIndex s2(5);

  s2.assign(s);
  BOOST_CHECK_EQUAL(s.getLevel(0), s2.getLevel(0));
  BOOST_CHECK_EQUAL(s.getIndex(1), s2.getIndex(1));
}

BOOST_AUTO_TEST_SUITE_END()
