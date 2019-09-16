// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifdef USE_SCALAPACK
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/datadriven/scalapack/BlacsProcessGrid.hpp>
#include <sgpp/datadriven/scalapack/DataMatrixDistributed.hpp>
#include <sgpp/datadriven/scalapack/DataVectorDistributed.hpp>

#include <mpi.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
using sgpp::datadriven::BlacsProcessGrid;
using sgpp::datadriven::DataMatrixDistributed;
using sgpp::datadriven::DataVectorDistributed;

struct FixtureDataVectorDistributed {
  FixtureDataVectorDistributed()
      : nrows(5),
        ncols(3),
        N(nrows * ncols),
        localGrid(std::make_shared<BlacsProcessGrid>(1, 1)),
        d_rand(localGrid, N, 2),
        d_rand_local(N),
        min(0),
        max(0),
        sum(0) {
    // setup vector for testing, same as in tests for DataMatrix / DataVector
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
        d_rand.set(i * ncols + j, l_rand[i][j]);
        d_rand_local.set(i * ncols + j, l_rand[i][j]);
      }
    }

    if (BlacsProcessGrid::availableProcesses() >= 2) {
      processGrid = std::make_shared<BlacsProcessGrid>(2, 1);
      d_rand_shared =
          std::make_shared<DataVectorDistributed>(d_rand_local.data(), processGrid, N, 2);
    }

    BOOST_TEST_MESSAGE("setup fixture");
  }
  ~FixtureDataVectorDistributed() {
    for (int i = 0; i < nrows; ++i) {
      delete[] l_rand[i];
    }

    delete[] l_rand;
    BOOST_TEST_MESSAGE("teardown fixture");
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
  std::shared_ptr<BlacsProcessGrid> localGrid;
  std::shared_ptr<BlacsProcessGrid> processGrid;
  std::shared_ptr<DataVectorDistributed> d_rand_shared;
  double** l_rand;
  DataVectorDistributed d_rand;
  DataVector d_rand_local;
  double min, max, sum;
};

BOOST_FIXTURE_TEST_SUITE(testDataVectorDistributed, FixtureDataVectorDistributed)

BOOST_AUTO_TEST_CASE(testConstructor) {
  int rowBlockSize = 10;
  DataVectorDistributed d(localGrid, 42, rowBlockSize);
  BOOST_CHECK_EQUAL(d.getGlobalRows(), 42);

  if (d.isProcessMapped()) {
    BOOST_CHECK_GT(d.getLocalRows(), 0);
    BOOST_CHECK_EQUAL(d.getMatrix().getLocalColumns(), 42);
    BOOST_CHECK_EQUAL(d.getMatrix().getLocalRows(), 1);

    for (size_t i = 0; i < d.getLocalRows(); i++) {
      BOOST_CHECK_CLOSE(d.getLocalPointer()[i], 0.0, 0.0001);
    }

  } else {
    BOOST_CHECK_EQUAL(d.getLocalRows(), 0);
    BOOST_CHECK_EQUAL(d.getMatrix().getLocalColumns(), 0);
  }

  int n = 9;
  int blockSize = 2;
  std::vector<double> testData{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
  DataVectorDistributed localVector(testData.data(), localGrid, n, blockSize);
  if (localGrid->isProcessInGrid()) {
    DataVector checkVector(testData.data(), n);

    // check that both matrices are equal
    for (int i = 0; i < n; i++) {
      BOOST_CHECK_CLOSE(checkVector.getPointer()[i], localVector.getLocalPointer()[i], 0.0001);
    }

    // test (local) set and get operations
    for (int i = 0; i < n; i++) {
      double value = i * 2;
      localVector.set(i, value);
      BOOST_CHECK_CLOSE(localVector.get(i), value, 0.0001);
      BOOST_CHECK_CLOSE(localVector.getLocalPointer()[i], value, 0.0001);
    }
  }

  if (processGrid) {
    std::vector<double> testData{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
    DataVectorDistributed testVector(testData.data(), processGrid, 9, 2);

    DataVector checkVector(testData.data(), 9);

    DataVector testVectorLocal = testVector.toLocalDataVector();

    // check that both matrices are equal
    if (testVector.isProcessMapped()) {
      for (int i = 0; i < 5; i++) {
        BOOST_CHECK_CLOSE(checkVector.get(i), testVector.get(i), 0.0001);
        if (processGrid->getCurrentProcess() == 0) {
          BOOST_CHECK_CLOSE(checkVector.get(i), testVectorLocal.get(i), 0.0001);
        }
      }
    }

    DataVectorDistributed vector(processGrid, 10, 2);
    // test distributed set and get operations
    for (int i = 0; i < 5; i++) {
      double initValue = vector.get(i);
      if (processGrid->isProcessInGrid()) {
        BOOST_CHECK_SMALL(initValue, 0.0001);
      }
      MPI_Barrier(MPI_COMM_WORLD);
      double value = i;
      vector.set(i, value);
      MPI_Barrier(MPI_COMM_WORLD);
      if (vector.isProcessMapped()) {
        double getValue = vector.get(i);
        BOOST_CHECK_CLOSE(getValue, value, 0.0001);
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(testOps) {
  double tolerance = 1e-12;

  DataVectorDistributed d = d_rand;
  double scalar = 5;

  // scale
  d = d_rand;
  if (d.isProcessMapped()) {
    d.scale(scalar);
    for (int i = 0; i < N; i++) {
      BOOST_CHECK_CLOSE(d(i), d_rand_local[i] * scalar, tolerance);
    }
  }

  if (processGrid && d_rand_shared->isProcessMapped()) {
    d = (*d_rand_shared);
    d.scale(scalar);
    for (int i = 0; i < N; i++) {
      BOOST_CHECK_CLOSE(d(i), d_rand_local[i] * scalar, tolerance);
    }
  }

  // dot
  d = d_rand;
  if (d.isProcessMapped()) {
    double dotTest = d.dot(d_rand);
    double dot = d_rand_local.dotProduct(d_rand_local);
    BOOST_CHECK_CLOSE(dotTest, dot, tolerance);
  }

  if (processGrid && d_rand_shared->isProcessMapped()) {
    d = (*d_rand_shared);
    double dotTest = d.dot(*d_rand_shared);
    double dot = d_rand_local.dotProduct(d_rand_local);
    BOOST_CHECK_CLOSE(dotTest, dot, tolerance);
  }

  // add
  d = d_rand;
  if (d.isProcessMapped()) {
    d.add(d_rand);
    for (int i = 0; i < N; i++) {
      BOOST_CHECK_CLOSE(d(i), d_rand_local[i] * 2, tolerance);
    }
  }

  if (processGrid && d_rand_shared->isProcessMapped()) {
    d = (*d_rand_shared);
    d.add(*d_rand_shared);
    for (int i = 0; i < N; i++) {
      BOOST_CHECK_CLOSE(d(i), d_rand_local[i] * 2, tolerance);
    }
  }
}

BOOST_AUTO_TEST_CASE(testSetUp) {
  if (localGrid->isProcessInGrid()) {
    for (int i = 0; i < nrows; ++i) {
      BOOST_CHECK_EQUAL(d_rand.get(i), d_rand_local.get(i));
    }
  }
}

BOOST_AUTO_TEST_CASE(testIndices) {
  int n = 12;
  int nb = 6;
  DataVectorDistributed local(localGrid, n, nb);

  if (local.isProcessMapped()) {
    for (int i = 0; i < n; i++) {
      BOOST_CHECK_EQUAL(local.localToGlobalRowIndex(i), i);
      BOOST_CHECK_EQUAL(local.globalToLocalRowIndex(i), i);
      BOOST_CHECK_EQUAL(local.globalToRowProcessIndex(i), 0);
    }
  }

  if (processGrid) {
    DataVectorDistributed test(processGrid, n, nb);
    int row = processGrid->getCurrentRow();

    if (test.isProcessMapped()) {
      for (int i = 0; i < n; i++) {
        if (i < nb) {
          BOOST_CHECK_EQUAL(test.localToGlobalRowIndex(i), row * nb + i);
        }
        BOOST_CHECK_EQUAL(test.globalToLocalRowIndex(i), i % nb);
        BOOST_CHECK_EQUAL(test.globalToRowProcessIndex(i), i / nb);
      }
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()

#endif
