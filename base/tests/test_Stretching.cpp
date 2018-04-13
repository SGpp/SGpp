// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <sgpp/base/grid/common/Stretching.hpp>

#include <string>
#include <vector>

using sgpp::base::BoundingBox1D;
using sgpp::base::DataVector;
using sgpp::base::Stretching;
using sgpp::base::Stretching1D;

BOOST_AUTO_TEST_SUITE(TestStretching)

BOOST_AUTO_TEST_CASE(testStretching) {
  const std::vector<BoundingBox1D> boundaries{{1.2, 3.4}, {-5.0, 5.0}};
  const std::vector<Stretching1D> stretching1Ds{Stretching1D("cc"), Stretching1D("id")};
  const Stretching stretching(boundaries, stretching1Ds);

  BOOST_CHECK_EQUAL(stretching.getStretching1D(0).type, "cc");
  BOOST_CHECK_EQUAL(stretching.getStretching1D(1).type, "id");
  BOOST_CHECK_EQUAL(stretching.getStretchingMode(), "analytic");

  BOOST_CHECK_EQUAL(stretching.getCoordinate(0, 0, 0), 1.2);
  BOOST_CHECK_EQUAL(stretching.getCoordinate(0, 1, 0), 3.4);
  BOOST_CHECK_CLOSE(stretching.getCoordinate(1, 1, 0), 2.3, 1e-8);
  BOOST_CHECK_CLOSE(stretching.getCoordinate(2, 1, 0), 1.522182540694798, 1e-8);
  BOOST_CHECK_CLOSE(stretching.getCoordinate(2, 3, 0), 3.077817459305202, 1e-8);
  BOOST_CHECK_CLOSE(stretching.getCoordinate(3, 1, 0), 1.283732514237585, 1e-8);
  BOOST_CHECK_CLOSE(stretching.getCoordinate(3, 3, 0), 1.879048224398401, 1e-8);
  BOOST_CHECK_CLOSE(stretching.getCoordinate(3, 5, 0), 2.720951775601599, 1e-8);
  BOOST_CHECK_CLOSE(stretching.getCoordinate(3, 7, 0), 3.316267485762416, 1e-8);

  BOOST_CHECK_EQUAL(stretching.getCoordinate(0, 0, 1), -5.0);
  BOOST_CHECK_EQUAL(stretching.getCoordinate(0, 1, 1), 5.0);
  BOOST_CHECK_EQUAL(stretching.getCoordinate(1, 1, 1), 0.0);
  BOOST_CHECK_EQUAL(stretching.getCoordinate(2, 1, 1), -2.5);
  BOOST_CHECK_EQUAL(stretching.getCoordinate(2, 3, 1), 2.5);
  BOOST_CHECK_EQUAL(stretching.getCoordinate(3, 1, 1), -3.75);
  BOOST_CHECK_EQUAL(stretching.getCoordinate(3, 3, 1), -1.25);
  BOOST_CHECK_EQUAL(stretching.getCoordinate(3, 5, 1), 1.25);
  BOOST_CHECK_EQUAL(stretching.getCoordinate(3, 7, 1), 3.75);

  const std::string serializedStretching = stretching.serialize();
  Stretching stretching2(2);
  stretching2.unserialize(serializedStretching, "analytic", SERIALIZATION_VERSION);
  BOOST_CHECK_EQUAL(stretching2.getIntervalOffset(0), 1.2);
  BOOST_CHECK_EQUAL(stretching2.getIntervalWidth(0), 2.2);
  BOOST_CHECK(!stretching2.hasDirichletBoundaryLeft(0));
  BOOST_CHECK(!stretching2.hasDirichletBoundaryRight(0));
  BOOST_CHECK_EQUAL(stretching2.getIntervalOffset(1), -5.0);
  BOOST_CHECK_EQUAL(stretching2.getIntervalWidth(1), 10.0);
  BOOST_CHECK(!stretching2.hasDirichletBoundaryLeft(1));
  BOOST_CHECK(!stretching2.hasDirichletBoundaryRight(1));
  BOOST_CHECK_EQUAL(stretching2.getStretching1D(0).type, "cc");
  BOOST_CHECK_EQUAL(stretching2.getStretching1D(1).type, "id");
  BOOST_CHECK_EQUAL(stretching2.getStretchingMode(), "analytic");
  BOOST_CHECK_CLOSE(stretching2.getCoordinate(3, 5, 0), 2.720951775601599, 1e-8);
  BOOST_CHECK_EQUAL(stretching2.getCoordinate(3, 5, 1), 1.25);
}

BOOST_AUTO_TEST_SUITE_END()
