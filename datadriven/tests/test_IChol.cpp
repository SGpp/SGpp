/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * test_IChol.cpp
 *
 *  Created on: Nov 26, 2016
 *      Author: Michael Lettrich
 */

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <sgpp/datadriven/algorithm/IChol.hpp>

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
using sgpp::datadriven::IChol;

BOOST_AUTO_TEST_SUITE(test_IChol)

BOOST_AUTO_TEST_CASE(decomp_identity) {
  auto size = 4u;
  DataMatrix A{size, size};

  // initialize
  A.setAll(0);
  for (auto i = 0u; i < size; i++) {
    A.set(i, i, 1);
  }

  auto B = A;

  // decomp:
  IChol::decompose(A, 1);

  // test
  for (auto i = 0u; i < A.getSize(); i++) {
    BOOST_CHECK_CLOSE(A[i], B[i], 10e-5);
  }
}

BOOST_AUTO_TEST_CASE(decomp_diag) {
  auto size = 4u;
  DataMatrix A{size, size};

  // initialize
  A.setAll(0);
  auto B = A;
  for (auto i = 0u; i < size; i++) {
    A.set(i, i, 4);
    B.set(i, i, 2);
  }

  // decomp:
  IChol::decompose(A, 1);

  // test
  for (auto i = 0u; i < A.getSize(); i++) {
    BOOST_CHECK_CLOSE(A[i], B[i], 10e-5);
  }
}

BOOST_AUTO_TEST_CASE(decomp_arbitrary) {
  auto size = 2u;

  // initialize
  double a_val[]{8.0, 4.0, 4.0, 9.0};
  double b_val[]{2.0 * sqrt(2), 4.0, sqrt(2), sqrt(7)};

  DataMatrix A{a_val, size, size};
  DataMatrix B{b_val, size, size};

  // decomp:
  IChol::decompose(A, 1);

  // test
  for (auto i = 0u; i < A.getSize(); i++) {
    BOOST_CHECK_CLOSE(A[i], B[i], 10e-5);
  }
}

BOOST_AUTO_TEST_CASE(norm) {
  auto size = 2u;

  // initialize
  double a_val[]{10.0, 5.0, 5.0, 20.0};
  double b_val[]{1.0, 5.0, 0.25, 1};
  double b_vec[]{10.0, 20.0};

  DataMatrix A{a_val, size, size};
  DataVector aNorm{size};
  DataMatrix B{b_val, size, size};
  DataVector bNorm{b_vec, size};

  // decomp:
  IChol::normToUnitDiagonal(A, aNorm);

  // test
  for (auto i = 0u; i < A.getSize(); i++) {
    BOOST_CHECK_CLOSE(A[i], B[i], 10e-5);
  }

  for (auto i = 0u; i < aNorm.getSize(); i++) {
    BOOST_CHECK_CLOSE(aNorm[i], bNorm[i], 10e-5);
  }
}

BOOST_AUTO_TEST_CASE(reaplyNorm) {
  auto size = 2u;

  // initialize
  double b_val[]{10.0, 5.0, 5.0, 20.0};
  double a_val[]{1.0, 5.0, 0.25, 1};
  double a_vec[]{10.0, 20.0};

  DataMatrix A{a_val, size, size};
  DataVector aNorm{a_vec, size};
  DataMatrix B{b_val, size, size};

  // decomp:
  IChol::reaplyDiagonal(A, aNorm);

  // test
  for (auto i = 0u; i < A.getSize(); i++) {
    BOOST_CHECK_CLOSE(A[i], B[i], 10e-5);
  }
}

BOOST_AUTO_TEST_SUITE_END()
