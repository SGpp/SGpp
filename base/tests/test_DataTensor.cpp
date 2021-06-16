// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK

#include <algorithm>
#include <cmath>
#include <iostream>

#include <boost/test/unit_test.hpp>

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataTensor.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

using sgpp::base::DataMatrix;
using sgpp::base::DataTensor;
using sgpp::base::DataVector;

BOOST_AUTO_TEST_SUITE(testDataTensor)

BOOST_AUTO_TEST_CASE(testTensorConstructor) {
  DataTensor t1(8, 42, 17);
  BOOST_CHECK_EQUAL(t1.getNdepth(), 8);
  BOOST_CHECK_EQUAL(t1.getNrows(), 42);
  BOOST_CHECK_EQUAL(t1.getNcols(), 17);

  DataTensor t2(8, 42, 17, 18.3);
  BOOST_CHECK_EQUAL(t2.getNdepth(), 8);
  BOOST_CHECK_EQUAL(t2.getNrows(), 42);
  BOOST_CHECK_EQUAL(t2.getNcols(), 17);
  BOOST_CHECK_EQUAL(t2.get(2, 3, 4), 18.3);
}

BOOST_AUTO_TEST_CASE(testGetSet) {
  DataTensor t1(4, 3, 2, 1.0);
  t1.set(2, 0, 1, 11.2);
  BOOST_CHECK_EQUAL(t1.get(2, 0, 0), 1.0);
  BOOST_CHECK_EQUAL(t1.get(2, 0, 1), 11.2);

  // std::cout << "t1:\n" << t1.toString() << "\n";

  DataMatrix m1;
  DataMatrix ref_m1(3, 2, 1.0);
  ref_m1.set(0, 1, 11.2);
  t1.getMatrix(2, m1);
  // std::cout << "m1:\n" << m1.toString() << "\n";
  // std::cout << "ref_m1:\n" << ref_m1.toString() << "\n";
  for (size_t i = 0; i < 3; i++) {
    for (size_t j = 0; j < 2; j++) {
      BOOST_CHECK_EQUAL(m1.get(i, j), ref_m1.get(i, j));
    }
  }

  DataVector v1;
  DataVector ref_v1(3, 1.0);
  ref_v1.set(0, 11.2);
  t1.getColumn(2, 1, v1);
  // std::cout << "v1:\n" << v1.toString() << "\n";
  // std::cout << "ref_v1:\n" << ref_v1.toString() << "\n";
  for (size_t i = 0; i < 3; i++) {
    BOOST_CHECK_EQUAL(v1[i], ref_v1[i]);
  }

  DataVector v2;
  DataVector ref_v2(2, 1.0);
  ref_v2.set(1, 11.2);
  t1.getRow(2, 0, v2);
  for (size_t j = 0; j < 2; j++) {
    BOOST_CHECK_EQUAL(v2[j], ref_v2[j]);
  }
}

BOOST_AUTO_TEST_SUITE_END()
