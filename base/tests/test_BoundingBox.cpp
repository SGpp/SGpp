// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <sgpp/base/grid/common/BoundingBox.hpp>

#include <string>

using sgpp::base::BoundingBox;
using sgpp::base::BoundingBox1D;
using sgpp::base::DataVector;

BOOST_AUTO_TEST_SUITE(TestBoundingBoxStretching)

BOOST_AUTO_TEST_CASE(testBoundingBox1D) {
  {
    const BoundingBox1D bb1D;
    BOOST_CHECK_EQUAL(bb1D.leftBoundary, 0.0);
    BOOST_CHECK_EQUAL(bb1D.rightBoundary, 1.0);
    BOOST_CHECK(!bb1D.bDirichletLeft);
    BOOST_CHECK(!bb1D.bDirichletRight);
  }

  {
    const BoundingBox1D bb1D(1.2, 3.4);
    BOOST_CHECK_EQUAL(bb1D.leftBoundary, 1.2);
    BOOST_CHECK_EQUAL(bb1D.rightBoundary, 3.4);
    BOOST_CHECK(!bb1D.bDirichletLeft);
    BOOST_CHECK(!bb1D.bDirichletRight);
  }

  {
    const BoundingBox1D bb1D(1.2, 3.4, true, true);
    BOOST_CHECK_EQUAL(bb1D.leftBoundary, 1.2);
    BOOST_CHECK_EQUAL(bb1D.rightBoundary, 3.4);
    BOOST_CHECK(bb1D.bDirichletLeft);
    BOOST_CHECK(bb1D.bDirichletRight);
  }
}

BOOST_AUTO_TEST_CASE(testBoundingBoxInOneDimension) {
  BoundingBox bb(3);
  BOOST_CHECK(bb.isUnitCube());

  const BoundingBox1D bb1D_1(1.2, 3.4);
  bb.setBoundary(1, bb1D_1);
  BOOST_CHECK(!bb.isUnitCube());
  BOOST_CHECK_EQUAL(bb.getIntervalOffset(1), 1.2);
  BOOST_CHECK_EQUAL(bb.getIntervalWidth(1), 2.2);
  BOOST_CHECK(!bb.hasDirichletBoundaryLeft(1));
  BOOST_CHECK(!bb.hasDirichletBoundaryRight(1));

  BOOST_CHECK_EQUAL(bb.transformPointToBoundingBox(1, 0.4), 2.08);
  BOOST_CHECK_EQUAL(bb.transformPointToUnitCube(1, 2.08), 0.4);
  BOOST_CHECK(bb.isContainingPoint(1, 2.08));
  BOOST_CHECK(bb.isContainingPoint(1, 1.2));
  BOOST_CHECK(bb.isContainingPoint(1, 3.4));
  BOOST_CHECK(!bb.isContainingPoint(1, 1.1));
  BOOST_CHECK(!bb.isContainingPoint(1, 3.5));

  const BoundingBox1D& bb1D_2 = bb.getBoundary(1);
  BOOST_CHECK_EQUAL(bb1D_1.leftBoundary, bb1D_2.leftBoundary);
  BOOST_CHECK_EQUAL(bb1D_1.rightBoundary, bb1D_2.rightBoundary);
  BOOST_CHECK_EQUAL(bb1D_1.bDirichletLeft, bb1D_2.bDirichletLeft);
  BOOST_CHECK_EQUAL(bb1D_1.bDirichletRight, bb1D_2.bDirichletRight);
}

BOOST_AUTO_TEST_CASE(testBoundingBoxInMultipleDimensions) {
  const BoundingBox bb({{1.2, 3.4, true, false}, {-5.0, 5.0, false, true}});
  BOOST_CHECK(!bb.isUnitCube());
  BOOST_CHECK_EQUAL(bb.getIntervalOffset(0), 1.2);
  BOOST_CHECK_EQUAL(bb.getIntervalWidth(0), 2.2);
  BOOST_CHECK(bb.hasDirichletBoundaryLeft(0));
  BOOST_CHECK(!bb.hasDirichletBoundaryRight(0));
  BOOST_CHECK_EQUAL(bb.getIntervalOffset(1), -5.0);
  BOOST_CHECK_EQUAL(bb.getIntervalWidth(1), 10.0);
  BOOST_CHECK(!bb.hasDirichletBoundaryLeft(1));
  BOOST_CHECK(bb.hasDirichletBoundaryRight(1));

  DataVector point(2);
  point[0] = 0.5;
  point[1] = 0.2;

  bb.transformPointToBoundingBox(point);
  BOOST_CHECK_CLOSE(point[0], 2.3, 1e-8);
  BOOST_CHECK_CLOSE(point[1], -3.0, 1e-8);
  BOOST_CHECK(bb.isContainingPoint(point));

  bb.transformPointToUnitCube(point);
  BOOST_CHECK_CLOSE(point[0], 0.5, 1e-8);
  BOOST_CHECK_CLOSE(point[1], 0.2, 1e-8);
  BOOST_CHECK(!bb.isContainingPoint(point));

  const std::string serializedBB = bb.serialize();
  BoundingBox bb2(2);
  bb2.unserialize(serializedBB, SERIALIZATION_VERSION);
  BOOST_CHECK_EQUAL(bb2.getIntervalOffset(0), 1.2);
  BOOST_CHECK_EQUAL(bb2.getIntervalWidth(0), 2.2);
  BOOST_CHECK(bb2.hasDirichletBoundaryLeft(0));
  BOOST_CHECK(!bb2.hasDirichletBoundaryRight(0));
  BOOST_CHECK_EQUAL(bb2.getIntervalOffset(1), -5.0);
  BOOST_CHECK_EQUAL(bb2.getIntervalWidth(1), 10.0);
  BOOST_CHECK(!bb2.hasDirichletBoundaryLeft(1));
  BOOST_CHECK(bb2.hasDirichletBoundaryRight(1));
}

BOOST_AUTO_TEST_SUITE_END()
