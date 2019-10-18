// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp/datadriven/algorithm/DBMatOfflineDenseIChol.hpp>

#include <omp.h>
#include <vector>

using sgpp::base::DataMatrix;

BOOST_AUTO_TEST_SUITE(DenseIChol_test)

BOOST_AUTO_TEST_CASE(decomp_identity) {
  // clang-format off
  const std::vector<double> data{
      1, 0, 0, 0, 0,
      0, 1, 0, 0, 0,
      0, 0, 1, 0, 0,
      0, 0, 0, 1, 0,
      0, 0, 0, 0, 1};
  // clang-format on

  auto size = 5u;
  DataMatrix A(data.data(), size, size);
  DataMatrix B = A;

  // decomp:
  sgpp::datadriven::DBMatOfflineDenseIChol::ichol(B, A, 1);

  // test
  for (auto i = 0u; i < A.getSize(); i++) {
    BOOST_CHECK_CLOSE(A[i], B[i], 1e-4);
  }
}

BOOST_AUTO_TEST_CASE(decomp_diag) {
  // clang-format off
  const std::vector<double> data{
      1, 0, 0, 0, 0,
      0, 4, 0, 0, 0,
      0, 0, 9, 0, 0,
      0, 0, 0, 16, 0,
      0, 0, 0, 0, 25};
  // clang-format on
  // clang-format off
  const std::vector<double> results{
      1, 0, 0, 0, 0,
      0, 2, 0, 0, 0,
      0, 0, 3, 0, 0,
      0, 0, 0, 4, 0,
      0, 0, 0, 0, 5};
  // clang-format on

  auto size = 5u;
  DataMatrix A(data, size);
  auto B = A;

  // decomp:
  sgpp::datadriven::DBMatOfflineDenseIChol::ichol(B, A, 1);

  // test
  for (auto i = 0u; i < A.getSize(); i++) {
    BOOST_CHECK_CLOSE(A[i], results[i], 1e-4);
  }
}

BOOST_AUTO_TEST_CASE(decomp_arbitrary) {
  // we only get reproducable results if we run on 1 omp thread
  auto numThreads = 0;

#pragma omp parallel
  {
#pragma omp single
    { numThreads = omp_get_num_threads(); }
  }
  omp_set_num_threads(1);

  auto size = 5u;

  // clang-format off
  const std::vector<double> data{
      2, 0, 0, 0, 0,
      2, 6, 0, 0, 0,
      1, 0, 6, 0, 0,
      0, 0, 1, 6, 0,
      1, 0, 4, 2, 16};
  // clang-format on
  // clang-format off
  const std::vector<double> results{
      1.414213562373095, 0, 0, 0, 0,
      1.414213562373095, 2.000000000000000, 0, 0, 0,
      0.707106781186547, 0, 2.345207879911715, 0, 0,
      0, 0, 0.426401432711221, 2.412090756622109 , 0,
      0.707106781186547, 0, 1.492405014489273, 0.565333771083307, 3.599045012221992};
  // clang-format on

  DataMatrix A(data, size);
  auto B = A;

  // decomp:
  sgpp::datadriven::DBMatOfflineDenseIChol::ichol(B, A, 1);

  // test
  for (auto i = 0u; i < A.getSize(); i++) {
    BOOST_CHECK_CLOSE(A[i], results[i], 1e-4);
  }
  omp_set_num_threads(numThreads);
}

BOOST_AUTO_TEST_SUITE_END()
