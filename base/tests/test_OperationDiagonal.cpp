// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/operation/hash/OperationDiagonal.hpp>

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
using sgpp::base::OperationDiagonal;

BOOST_AUTO_TEST_SUITE(TestOperationDiagonal)

BOOST_AUTO_TEST_CASE(testOperationDiagonal) {
  /*
  double testArr[] = {0.0, 1.0, 2.0};

  auto identVec = DataVector(3, 1.0);
  auto testVec = DataVector(testArr, 3);
  auto zeroVec = DataVector(3, 0.0);

  auto identOp = OperationDiagonal(identVec);
  auto testOp = OperationDiagonal(testVec);
  auto zeroOp = OperationDiagonal(zeroVec);

  auto is1 = DataVector(3);
  zeroOp.mult(testVec, is1);
  auto should1 = 0.0;
  BOOST_CHECK_CLOSE(is1.l2Norm(), should1, 1e-7);

  auto is2 = DataVector(3);
  identOp.mult(testVec, is2);
  auto should2 = testVec;
  BOOST_CHECK_CLOSE(is2.l2Norm(), should2.l2Norm(), 1e-7);

  auto is3 = DataVector(3);
  testOp.mult(testVec, is3);
  double should3Arr[] = {0.0, 1.0, 4.0};
  auto should3 = DataVector(should3Arr, 3);
  BOOST_CHECK_CLOSE(is3.l2Norm(), should3.l2Norm(), 1e-7);
  */
  BOOST_CHECK_CLOSE(9000.0, 42.0, 0.0);  // TODO(krenzls) write test!
}

BOOST_AUTO_TEST_SUITE_END()
