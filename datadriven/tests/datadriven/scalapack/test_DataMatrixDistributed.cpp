/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * test_DataMatrixDistributed.cpp
 *
 * Created on: Feb 19, 2019
 *     Author: Jan Schopohl
 */

#ifdef USE_SCALAPACK
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <sgpp/datadriven/scalapack/BlacsProcessGrid.hpp>
#include <sgpp/datadriven/scalapack/DataMatrixDistributed.hpp>
#include <sgpp/datadriven/scalapack/DataVectorDistributed.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>

using sgpp::datadriven::BlacsProcessGrid;
using sgpp::datadriven::DataMatrixDistributed;
using sgpp::datadriven::DataVectorDistributed;

struct FixtureDataMatrixDistributed {
  FixtureDataMatrixDistributed()
      : nrows(5),
        ncols(3),
        N(nrows * ncols),
        processGrid(std::make_shared<BlacsProcessGrid>(1, 1)),
        d_rand(processGrid, nrows, ncols, DataMatrixDistributed::DTYPE::DENSE, 2, 2),
        min(0),
        max(0),
        sum(0) {
    // TODO(jan) size of the process grid?

    // setup matrix for testing, same as in tests for DataMatrix
    l_rand = new double*[nrows];

    for (int i = 0; i < nrows; ++i) {
      l_rand[i] = new double[ncols];
    }

    for (int i = 0; i < nrows; ++i) {
      for (int j = 0; j < ncols; ++j) {
        l_rand[i][j] = i * j + i * 0.5 + 2.34 * j;
        min = min > l_rand[i][j] ? l_rand[i][j] : min;
        max = max < l_rand[i][j] ? l_rand[i][j] : max;
        sum += l_rand[i][j];
      }
    }

    for (int i = 0; i < nrows; ++i) {
      for (int j = 0; j < ncols; ++j) {
        d_rand.set(i, j, l_rand[i][j]);
      }
    }

    BOOST_TEST_MESSAGE("setup fixture");
  }
  ~FixtureDataMatrixDistributed() {
    for (int i = 0; i < nrows; ++i) {
      delete[] l_rand[i];
    }

    delete[] l_rand;
    BOOST_TEST_MESSAGE("teardown fixture");
  }
  int nrows, ncols, N;
  std::shared_ptr<BlacsProcessGrid> processGrid;
  double** l_rand;
  DataMatrixDistributed d_rand;
  double min, max, sum;
};

BOOST_FIXTURE_TEST_SUITE(testDataMatrixDistributed, FixtureDataMatrixDistributed)

BOOST_AUTO_TEST_CASE(testConstructor) {
  DataMatrix d(42, 17);
  BOOST_CHECK_EQUAL(d.getSize(), 42 * 17);
  BOOST_CHECK_EQUAL(d.getNrows(), 42);
  BOOST_CHECK_EQUAL(d.getNcols(), 17);

  int nrows, ncols = 3;
  std::vector<double> testData{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
  // std::shared_ptr<BlacsProcessGrid> localGrid = std::make_shared<BlacsProcessGrid>(1, 1);
  DataMatrixDistributed localMatrix(testData.data(), processGrid, nrows, ncols,
                                    DataMatrixDistributed::DTYPE::DENSE, 2, 2);

  DataMatrix checkMatrix(testData.data(), 3, 3);

  // check that both matrices are equal
  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < ncols; j++) {
      BOOST_CHECK_EQUAL(checkMatrix.getPointer()[(i * ncols) + j],
                        localMatrix.getLocalPointer()[(i * ncols) + j]);
    }
  }

  // test if matrix is correctly created, test set and get operations
}

BOOST_AUTO_TEST_CASE(testSetUp) {
  for (int i = 0; i < nrows; ++i) {
    for (int j = 0; j < ncols; ++j) {
      BOOST_CHECK_EQUAL(d_rand.get(i, j), l_rand[i][j]);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()

#endif
