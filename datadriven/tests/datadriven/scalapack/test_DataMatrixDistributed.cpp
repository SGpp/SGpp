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

#include <mpi.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/datadriven/scalapack/BlacsProcessGrid.hpp>
#include <sgpp/datadriven/scalapack/DataMatrixDistributed.hpp>
#include <sgpp/datadriven/scalapack/DataVectorDistributed.hpp>

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
using sgpp::datadriven::BlacsProcessGrid;
using sgpp::datadriven::DataMatrixDistributed;
using sgpp::datadriven::DataVectorDistributed;

struct FixtureBlacsGrid {
  FixtureBlacsGrid() { BlacsProcessGrid::initializeBlacs(); }

  ~FixtureBlacsGrid() { BlacsProcessGrid::exitBlacs(); }
};

BOOST_GLOBAL_FIXTURE(FixtureBlacsGrid);

struct FixtureDataMatrixDistributed {
  FixtureDataMatrixDistributed()
      : nrows(5),
        ncols(3),
        N(nrows * ncols),
        processGrid(std::make_shared<BlacsProcessGrid>(1, 1)),
        d_rand(processGrid, nrows, ncols, 2, 2),
        d_rand_local(nrows, ncols),
        min(0),
        max(0),
        sum(0) {
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
        d_rand_local.set(i, j, l_rand[i][j]);
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

  static void assertMatrixClose(const DataMatrix& expected, const DataMatrixDistributed& actual) {
    if (actual.isProcessMapped()) {
      BOOST_CHECK_EQUAL(expected.getNrows(), actual.getGlobalRows());
      BOOST_CHECK_EQUAL(expected.getNcols(), actual.getGlobalCols());

      for (size_t i = 0; i < actual.getGlobalRows(); i++) {
        for (size_t j = 0; j < actual.getGlobalCols(); j++) {
          BOOST_CHECK_CLOSE(expected(i, j), actual(i, j), 0.0001);
        }
      }
    }
  }

  static void assertVectorClose(const DataVector& expected, const DataVectorDistributed& actual) {
    if (actual.isProcessMapped()) {
      BOOST_CHECK_EQUAL(expected.getSize(), actual.getGlobalRows());

      for (size_t i = 0; i < actual.getGlobalRows(); i++) {
        BOOST_CHECK_CLOSE(expected.get(i), actual(i), 0.0001);
      }
    }
  }

  int nrows, ncols, N;
  std::shared_ptr<BlacsProcessGrid> processGrid;
  double** l_rand;
  DataMatrixDistributed d_rand;
  DataMatrix d_rand_local;
  double min, max, sum;
};

BOOST_FIXTURE_TEST_SUITE(testDataMatrixDistributed, FixtureDataMatrixDistributed)

BOOST_AUTO_TEST_CASE(testConstructor) {
  int rowBlockSize = 10;
  int columnBlockSize = 1;
  DataMatrixDistributed d(processGrid, 42, 17, columnBlockSize, rowBlockSize);
  BOOST_CHECK_EQUAL(d.getGlobalRows(), 42);
  BOOST_CHECK_EQUAL(d.getGlobalCols(), 17);

  if (d.isProcessMapped()) {
    BOOST_CHECK_GT(d.getLocalRows(), 0);
    BOOST_CHECK_GT(d.getLocalRows(), 0);

    for (int i = 0; i < d.getLocalRows(); i++) {
      for (int j = 0; j < d.getLocalColumns(); j++) {
        BOOST_CHECK_CLOSE(d.getLocalPointer()[(i * 17) + j], 0.0, 0.0001);
      }
    }

  } else {
    BOOST_CHECK_EQUAL(d.getLocalRows(), 0);
    BOOST_CHECK_EQUAL(d.getLocalColumns(), 0);
  }

  std::shared_ptr<BlacsProcessGrid> localGrid = std::make_shared<BlacsProcessGrid>(1, 1);

  int n = 3;
  int blockSize = 2;
  std::vector<double> testData{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
  DataMatrixDistributed localMatrix(testData.data(), localGrid, n, n, blockSize, blockSize);
  if (localGrid->isProcessInGrid()) {
    DataMatrix checkMatrix(testData.data(), n, n);

    // check that both matrices are equal
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        BOOST_CHECK_CLOSE(checkMatrix.getPointer()[(i * n) + j],
                          localMatrix.getLocalPointer()[(i * n) + j], 0.0001);
      }
    }

    // test (local) set and get operations
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        double value = (i + 2) * (j + 2);
        localMatrix.set(i, j, value);
        BOOST_CHECK_CLOSE(localMatrix.get(i, j), value, 0.0001);
        BOOST_CHECK_CLOSE(localMatrix.getLocalPointer()[(i * n) + j], value, 0.0001);
      }
    }
  }

  if (BlacsProcessGrid::availableProcesses() >= 4) {
    std::shared_ptr<BlacsProcessGrid> grid = std::make_shared<BlacsProcessGrid>(2, 2);

    std::vector<double> testData{1.0,  2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,  9.0,
                                 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0,
                                 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0};
    DataMatrixDistributed testMatrix(testData.data(), grid, 5, 5, 2, 2);

    DataMatrix checkMatrix(testData.data(), 5, 5);

    DataMatrix testMatrixLocal = testMatrix.toLocalDataMatrix();

    // check that both matrices are equal
    if (testMatrix.isProcessMapped()) {
      for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
          BOOST_CHECK_CLOSE(checkMatrix.get(i, j), testMatrix.get(i, j), 0.0001);
          if (grid->getCurrentProcess() == 0) {
            BOOST_CHECK_CLOSE(checkMatrix.get(i, j), testMatrixLocal.get(i, j), 0.0001);
          }
        }
      }
    }

    DataMatrixDistributed matrix(grid, 5, 5, 2, 2);
    // test distributed set and get operations
    for (int i = 0; i < 5; i++) {
      for (int j = 0; j < 5; j++) {
        double initValue = matrix.get(i, j);
        if (grid->isProcessInGrid()) {
          BOOST_CHECK_SMALL(initValue, 0.0001);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        double value = (i + 1) * (j + 1);
        matrix.set(i, j, value);
        MPI_Barrier(MPI_COMM_WORLD);
        if (matrix.isProcessMapped()) {
          double getValue = matrix.get(i, j);
          BOOST_CHECK_CLOSE(getValue, value, 0.0001);
        }
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(testSetUp) {
  if (processGrid->isProcessInGrid()) {
    for (int i = 0; i < nrows; ++i) {
      for (int j = 0; j < ncols; ++j) {
        BOOST_CHECK_EQUAL(d_rand.get(i, j), l_rand[i][j]);
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(testTranspose) {
  DataMatrixDistributed d_rand_transposed = d_rand.transpose();
  d_rand_local.transpose();

  assertMatrixClose(d_rand_local, d_rand_transposed);

  if (processGrid->availableProcesses() >= 4) {
    std::shared_ptr<BlacsProcessGrid> grid = std::make_shared<BlacsProcessGrid>(2, 2);
    DataMatrixDistributed a(d_rand_local.data(), grid, ncols, nrows, 2, 2);

    DataMatrixDistributed a_transposed = a.transpose();
    d_rand_local.transpose();

    assertMatrixClose(d_rand_local, a_transposed);
  }
}

BOOST_AUTO_TEST_CASE(testPdgemv) {
  std::vector<double> testVectorData;
  for (int i = 0; i < ncols; i++) {
    testVectorData.push_back(i * 1.25);
  }

  DataVectorDistributed x(testVectorData.data(), processGrid, ncols, 2);

  DataVectorDistributed result(processGrid, nrows, 2);

  DataMatrixDistributed::mult(d_rand, x, result);

  DataVector expectedX(testVectorData.data(), ncols);
  DataVector expected(nrows);

  d_rand_local.mult(expectedX, expected);
  assertVectorClose(expected, result);

  std::vector<double> testData(nrows, 1.0);
  result = DataVectorDistributed(testData.data(), processGrid, nrows, 2);

  DataMatrixDistributed::mult(d_rand, x, result, false, 1.0, 1.0);

  DataVector offset(testData.data(), nrows);
  expected.add(offset);

  assertVectorClose(expected, result);

  if (BlacsProcessGrid::availableProcesses() >= 4) {
    std::shared_ptr<BlacsProcessGrid> grid = std::make_shared<BlacsProcessGrid>(2, 2);
    DataMatrixDistributed a(d_rand_local.data(), grid, nrows, ncols, 2, 2);
    DataVectorDistributed x(testVectorData.data(), grid, ncols, 2);
    DataVectorDistributed result(grid, nrows, 2);

    DataMatrixDistributed::mult(a, x, result);

    assertVectorClose(expected, result);
  }
}

BOOST_AUTO_TEST_CASE(testPdgemm) {
  std::vector<double> testVecA{1.0, 0.0};
  std::vector<double> testVecB{1.0, 0.0};
  std::vector<double> resultVecData{1.0};
  std::vector<double> resultVecDataT{1.0, 0.0, 0.0, 0.0};

  // test with vectors, note that different block sizes for rows and columns are also tested
  DataMatrixDistributed av(testVecA.data(), processGrid, 1, 2, 1, 2);
  DataMatrixDistributed bv(testVecB.data(), processGrid, 2, 1, 1, 2);
  DataMatrixDistributed cv(processGrid, 1, 1, 1, 2);
  DataMatrix resultVec(resultVecData.data(), 1, 1);

  DataMatrixDistributed::mult(av, bv, cv);

  assertMatrixClose(resultVec, cv);

  DataMatrixDistributed avt(testVecA.data(), processGrid, 2, 1, 2, 2);
  DataMatrixDistributed bvt(testVecB.data(), processGrid, 1, 2, 2, 2);
  DataMatrixDistributed cvt(processGrid, 2, 2, 2, 2);
  DataMatrix resultVecT(resultVecDataT.data(), 2, 2);

  DataMatrixDistributed::mult(avt, bvt, cvt);

  assertMatrixClose(resultVecT, cvt);

  std::vector<double> testA{1.0, 0.0, 1.0, 0.0};
  std::vector<double> testB{1.0, 0.0, 1.0, 1.0};
  std::vector<double> resultData{1.0, 0.0, 1.0, 0.0};

  DataMatrixDistributed a(testA.data(), processGrid, 2, 2, 2, 2);
  DataMatrixDistributed b(testB.data(), processGrid, 2, 2, 2, 2);

  DataMatrix resultCheck(resultData.data(), 2, 2);

  DataMatrixDistributed c(processGrid, 2, 2, 2, 2);

  DataMatrixDistributed::mult(a, b, c);

  assertMatrixClose(resultCheck, c);

  if (BlacsProcessGrid::availableProcesses() >= 4) {
    std::vector<double> testDataA{1.0, 0.0, 1.0, 0.0, 0.0, 2.0, 0.0, 1.0, 1.0, 3.0, 0.0, 0.0};
    std::vector<double> testDataB{1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 2.0, 2.0, 0.0, 0.0, 1.0};

    std::vector<double> resultData2{1.0, 2.0, 3.0, 2.0, 2.0, 3.0, 4.0, 3.0, 4.0};

    std::shared_ptr<BlacsProcessGrid> grid = std::make_shared<BlacsProcessGrid>(2, 2);

    std::vector<double> testASquare{1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0};
    std::vector<double> testBSquare{1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    DataMatrixDistributed as(testASquare.data(), grid, 3, 3, 2, 2);
    DataMatrixDistributed bs(testBSquare.data(), grid, 3, 3, 2, 2);
    DataMatrixDistributed cs(grid, 3, 3, 2, 2);

    DataMatrixDistributed::mult(as, bs, cs);

    DataMatrixDistributed a(testDataA.data(), grid, 3, 4, 2, 2);
    DataMatrixDistributed b(testDataB.data(), grid, 4, 3, 2, 2);
    DataMatrixDistributed c(grid, 3, 3, 2, 2);
    DataMatrix result(resultData2.data(), 3, 3);

    DataMatrixDistributed::mult(a, b, c);

    assertMatrixClose(result, c);
  }
}

BOOST_AUTO_TEST_CASE(testPgdeadd) {
  std::vector<double> testData;
  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < ncols; j++) {
      testData.push_back((i + j) * 1.25);
    }
  }

  DataMatrixDistributed x(testData.data(), processGrid, nrows, ncols, 2, 2);

  DataMatrixDistributed::add(d_rand, x);

  DataMatrix localX(testData.data(), nrows, ncols);

  d_rand_local.add(localX);

  assertMatrixClose(d_rand_local, d_rand);

  d_rand.sub(x);
  d_rand_local.sub(localX);

  assertMatrixClose(d_rand_local, d_rand);

  if (BlacsProcessGrid::availableProcesses() >= 4) {
    std::shared_ptr<BlacsProcessGrid> grid = std::make_shared<BlacsProcessGrid>(2, 2);
    DataMatrixDistributed c(d_rand_local.data(), grid, nrows, ncols, 2, 2);
    DataMatrixDistributed x(testData.data(), grid, nrows, ncols, 2, 2);

    DataMatrixDistributed::add(c, x);

    d_rand_local.add(localX);

    assertMatrixClose(d_rand_local, c);

    c.sub(x);
    d_rand_local.sub(localX);

    assertMatrixClose(d_rand_local, c);
  }
}

BOOST_AUTO_TEST_SUITE_END()

#endif
