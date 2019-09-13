// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp/base/datatypes/DataVector.hpp>

#include <algorithm>
#include <cmath>

using sgpp::base::DataVector;

struct FixtureDataVector {
  FixtureDataVector() : nrows(5), ncols(3), N(nrows * ncols), d_rand(N), min(0), max(0), sum(0) {
    l_rand_total = new double[nrows * ncols];
    l_rand = new double* [nrows];

    for (int i = 0; i < nrows; ++i) {
      l_rand[i] = new double[ncols];
    }

    for (int i = 0; i < nrows; ++i) {
      for (int j = 0; j < ncols; ++j) {
        l_rand[i][j] = i * j + i * 0.5 + 2.34 * j;
        l_rand_total[i * ncols + j] = l_rand[i][j];
        min = min > l_rand_total[i * ncols + j] ? l_rand_total[i * ncols + j] : min;
        max = max < l_rand_total[i * ncols + j] ? l_rand_total[i * ncols + j] : max;
        sum += l_rand_total[i * ncols + j];
      }
    }

    for (int i = 0; i < N; ++i) {
      d_rand[i] = l_rand_total[i];
    }

    BOOST_TEST_MESSAGE("setup fixture");
  }
  ~FixtureDataVector() {
    delete[] l_rand_total;

    for (int i = 0; i < nrows; ++i) {
      delete[] l_rand[i];
    }

    delete[] l_rand;
    BOOST_TEST_MESSAGE("teardown fixture");
  }
  int nrows, ncols, N;
  double** l_rand;
  double* l_rand_total;
  DataVector d_rand;
  double min, max, sum;
};

BOOST_FIXTURE_TEST_SUITE(testDataVector, FixtureDataVector)

BOOST_AUTO_TEST_CASE(testConstructor) {
  DataVector d = DataVector(2);
  BOOST_CHECK_EQUAL(d.getSize(), 2U);

  d = DataVector({1.0, 0.5, -2.0});
  BOOST_CHECK_EQUAL(d.getSize(), 3U);
  BOOST_CHECK_EQUAL(d[0], 1.0);
  BOOST_CHECK_EQUAL(d[1], 0.5);
  BOOST_CHECK_EQUAL(d[2], -2.0);
}

BOOST_AUTO_TEST_CASE(testSetUp) {
  for (int i = 0; i < N; ++i) {
    BOOST_CHECK_EQUAL(l_rand_total[i], d_rand[i]);
  }
}

BOOST_AUTO_TEST_CASE(testMinMax) {
  BOOST_CHECK_EQUAL(d_rand.min(), min);
  BOOST_CHECK_EQUAL(d_rand.max(), max);
}

BOOST_AUTO_TEST_CASE(testOps) {
  double tol = 1e-12;

  DataVector d = d_rand;
  DataVector d2(N);
  double scalar = 0.213;

  for (int i = 0; i < N; ++i) {
    d2[i] = static_cast<double>(i);
  }

  // add
  d = DataVector(d_rand);
  d.add(d2);
  for (int i = 0; i < N; ++i) {
    BOOST_CHECK_EQUAL(d[i], d_rand[i] + static_cast<double>(i));
  }

  // axpy
  d = DataVector(d_rand);
  d.axpy(scalar, d2);
  for (int i = 0; i < N; ++i) {
    BOOST_CHECK_EQUAL(d[i], d_rand[i] + scalar * static_cast<double>(i));
  }

  // dotProduct
  d = DataVector(d_rand);
  double dotProdResultActual = d.dotProduct(d2);
  double dotProdResultExact = 0;
  for (int i = 0; i < N; ++i) {
    dotProdResultExact += d_rand[i] * static_cast<double>(i);
  }
  BOOST_CHECK_EQUAL(dotProdResultActual, dotProdResultExact);

  // L2Norm
  d = DataVector(d_rand);
  double lTwoNormSquaredActual = d.l2Norm() * d.l2Norm();
  double lTwoNormSquaredExact = 0;
  for (int i = 0; i < N; ++i) {
    lTwoNormSquaredExact += d_rand[i] * d_rand[i];
  }
  BOOST_CHECK_CLOSE(lTwoNormSquaredActual, lTwoNormSquaredExact, tol);

  // max
  d = DataVector(d_rand);
  double maxActual = d.max();
  BOOST_CHECK_EQUAL(max, maxActual);

  // maxNorm
  d = DataVector(d_rand);
  double maxNormActual = d.maxNorm();
  double maxNormExpected = 0.0;
  for (int i = 1; i < N; ++i) {
    maxNormExpected = maxNormExpected < std::abs(d_rand[i]) ? std::abs(d_rand[i]) : maxNormExpected;
  }
  BOOST_CHECK_CLOSE(maxNormActual, maxNormExpected, tol);

  // min
  d = DataVector(d_rand);
  double minActual = d.min();
  BOOST_CHECK_EQUAL(minActual, min);

  // minmax
  d.minmax(&minActual, &maxActual);
  BOOST_CHECK_EQUAL(minActual, min);
  BOOST_CHECK_EQUAL(maxActual, max);

  // normalize
  d = DataVector(d_rand);
  d.normalize();
  double border = 0.0;
  double delta = (d_rand.max() - d_rand.min()) / (1 - 2 * border);
  for (int i = 0; i < N; i++) {
    BOOST_CHECK_CLOSE(d[i], (d_rand[i] - d_rand.min()) / delta + border, tol);
  }

  // normalize with border
  d = DataVector(d_rand);
  border = 3.64;
  d.normalize(border);
  delta = (d_rand.max() - d_rand.min()) / (1 - 2 * border);
  for (int i = 0; i < N; i++) {
    BOOST_CHECK_CLOSE(d[i], (d_rand[i] - d_rand.min()) / delta + border, tol);
  }

  // RMSNorm
  d = DataVector(d_rand);
  double rmsNormSquaredActual = d.RMSNorm() * d.RMSNorm();
  double rmsNormSquaredExpected = 0;
  for (int i = 0; i < N; ++i) {
    rmsNormSquaredExpected += d_rand[i] * d_rand[i];
  }
  rmsNormSquaredExpected *= 1.0 / N;
  BOOST_CHECK_CLOSE(rmsNormSquaredActual, rmsNormSquaredExpected, tol);

  // sub
  d = DataVector(d_rand);
  d.sub(d2);
  for (int i = 0; i < N; ++i) {
    BOOST_CHECK_EQUAL(d[i], d_rand[i] - static_cast<double>(i));
  }

  // mult scalar
  d = DataVector(d_rand);
  d.mult(scalar);
  for (int i = 0; i < N; ++i) {
    BOOST_CHECK_EQUAL(d[i], d_rand[i] * scalar);
  }

  // sum
  d = DataVector(d_rand);
  BOOST_CHECK_CLOSE(d.sum(), sum, tol);

  // square
  d = DataVector(d_rand);
  d.sqr();

  for (int i = 0; i < N; ++i) {
    BOOST_CHECK_EQUAL(d[i], d_rand[i] * d_rand[i]);
  }

  // sqrt
  d = DataVector(d_rand);
  d.sqr();
  d.sqrt();
  for (int i = 0; i < N; ++i) {
    BOOST_CHECK_CLOSE(d[i], d_rand[i], tol);
  }

  // abs
  d = DataVector(d_rand);
  d.abs();

  for (int i = 0; i < N; ++i) {
    BOOST_CHECK_EQUAL(d[i], std::abs(d_rand[i]));
  }

  // componentwise mult
  d = DataVector(d_rand);

  d.componentwise_mult(d2);

  for (int i = 0; i < N; ++i) {
    BOOST_CHECK_EQUAL(d[i], d_rand[i] * static_cast<double>(i));
  }

  // componentwise div
  d = DataVector(d_rand);

  for (int i = 0; i < N; ++i) {
    d2[i] = i + 1.0;
  }

  d.componentwise_div(d2);

  for (int i = 0; i < N; ++i) {
    BOOST_CHECK_CLOSE(d_rand[i] / (i + 1.0), d[i], 1e-12);
  }
}

BOOST_AUTO_TEST_CASE(testDotProduct) {
  DataVector d = DataVector(3);
  double x = 0;

  for (unsigned int i = 0; i < d.getSize(); ++i) {
    d[i] = static_cast<double>(i + 1);
    x += d[i] * d[i];
  }

  BOOST_CHECK_EQUAL(d.dotProduct(d), x);
}

BOOST_AUTO_TEST_SUITE_END()
