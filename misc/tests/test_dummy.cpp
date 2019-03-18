// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(TestDummy)

BOOST_AUTO_TEST_CASE(testDummy1) {
  BOOST_CHECK_EQUAL(1, 1);
}

BOOST_AUTO_TEST_SUITE_END()
