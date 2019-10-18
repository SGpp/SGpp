// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp/base/datatypes/DataVectorSP.hpp>

#include <algorithm>
#include <cmath>

using sgpp::base::DataVectorSP;

struct FixtureDataVectorSP {
  FixtureDataVectorSP() : nrows(5), ncols(3), N(nrows * ncols), d_rand(N), min(0), max(0), sum(0) {
    l_rand_total = new float[nrows * ncols];
    l_rand = new float*[nrows];

    for (int i = 0; i < nrows; ++i) {
      l_rand[i] = new float[ncols];
    }

    for (int i = 0; i < nrows; ++i) {
      for (int j = 0; j < ncols; ++j) {
        l_rand[i][j] = static_cast<float>(i) * static_cast<float>(j) +
                       static_cast<float>(i) * 0.5f + 2.34f * static_cast<float>(j);
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
  ~FixtureDataVectorSP() {
    delete[] l_rand_total;

    for (int i = 0; i < nrows; ++i) {
      delete[] l_rand[i];
    }

    delete[] l_rand;
    BOOST_TEST_MESSAGE("teardown fixture");
  }
  int nrows, ncols, N;
  float** l_rand;
  float* l_rand_total;
  DataVectorSP d_rand;
  float min, max, sum;
};

BOOST_FIXTURE_TEST_SUITE(testDataVectorSP, FixtureDataVectorSP)

BOOST_AUTO_TEST_CASE(testConstructor) {
  DataVectorSP d = DataVectorSP(2);
  BOOST_CHECK_EQUAL(d.getSize(), 2U);

  d = DataVectorSP({1.0f, 0.5f, -2.0f});
  BOOST_CHECK_EQUAL(d.getSize(), 3U);
  BOOST_CHECK_EQUAL(d[0], 1.0f);
  BOOST_CHECK_EQUAL(d[1], 0.5f);
  BOOST_CHECK_EQUAL(d[2], -2.0f);
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
  float tol = 1e-4f;

  DataVectorSP d = d_rand;
  DataVectorSP d2(N);
  float scalar = 0.213f;

  for (int i = 0; i < N; ++i) {
    d2[i] = static_cast<float>(i);
  }

  // add
  d = DataVectorSP(d_rand);
  d.add(d2);
  for (int i = 0; i < N; ++i) {
    BOOST_CHECK_EQUAL(d[i], d_rand[i] + static_cast<float>(i));
  }

  // axpy
  d = DataVectorSP(d_rand);
  d.axpy(scalar, d2);
  for (int i = 0; i < N; ++i) {
    BOOST_CHECK_EQUAL(d[i], d_rand[i] + scalar * static_cast<float>(i));
  }

  // dotProduct
  d = DataVectorSP(d_rand);
  float dotProdResultActual = d.dotProduct(d2);
  float dotProdResultExact = 0.0f;
  for (int i = 0; i < N; ++i) {
    dotProdResultExact += d_rand[i] * static_cast<float>(i);
  }
  BOOST_CHECK_EQUAL(dotProdResultActual, dotProdResultExact);

  // L2Norm
  d = DataVectorSP(d_rand);
  float lTwoNormSquaredActual = d.l2Norm() * d.l2Norm();
  float lTwoNormSquaredExact = 0.0f;
  for (int i = 0; i < N; ++i) {
    lTwoNormSquaredExact += d_rand[i] * d_rand[i];
  }
  BOOST_CHECK_CLOSE(lTwoNormSquaredActual, lTwoNormSquaredExact, tol);

  // max
  d = DataVectorSP(d_rand);
  float maxActual = d.max();
  BOOST_CHECK_EQUAL(max, maxActual);

  // maxNorm
  d = DataVectorSP(d_rand);
  float maxNormActual = d.maxNorm();
  float maxNormExpected = 0.0f;
  for (int i = 1; i < N; ++i) {
    maxNormExpected = maxNormExpected < static_cast<float>(std::abs(d_rand[i]))
                          ? static_cast<float>(std::abs(d_rand[i]))
                          : maxNormExpected;
  }
  BOOST_CHECK_CLOSE(maxNormActual, maxNormExpected, tol);

  // min
  d = DataVectorSP(d_rand);
  float minActual = d.min();
  BOOST_CHECK_EQUAL(minActual, min);

  // minmax
  d.minmax(&minActual, &maxActual);
  BOOST_CHECK_EQUAL(minActual, min);
  BOOST_CHECK_EQUAL(maxActual, max);

  // normalize
  d = DataVectorSP(d_rand);
  d.normalize();
  float border = 0.0f;
  float delta = (d_rand.max() - d_rand.min()) / (1.0f - 2.0f * border);
  for (int i = 0; i < N; i++) {
    BOOST_CHECK_CLOSE(d[i], (d_rand[i] - d_rand.min()) / delta + border, tol);
  }

  // normalize with border
  d = DataVectorSP(d_rand);
  border = 3.64f;
  d.normalize(border);
  delta = (d_rand.max() - d_rand.min()) / (1.0f - 2.0f * border);
  for (int i = 0; i < N; i++) {
    BOOST_CHECK_CLOSE(d[i], (d_rand[i] - d_rand.min()) / delta + border, tol);
  }

  // RMSNorm
  d = DataVectorSP(d_rand);
  float rmsNormSquaredActual = d.RMSNorm() * d.RMSNorm();
  float rmsNormSquaredExpected = 0.0f;
  for (int i = 0; i < N; ++i) {
    rmsNormSquaredExpected += d_rand[i] * d_rand[i];
  }
  rmsNormSquaredExpected *= 1.0f / static_cast<float>(N);
  BOOST_CHECK_CLOSE(rmsNormSquaredActual, rmsNormSquaredExpected, tol);

  // sub
  d = DataVectorSP(d_rand);
  d.sub(d2);
  for (int i = 0; i < N; ++i) {
    BOOST_CHECK_EQUAL(d[i], d_rand[i] - static_cast<float>(i));
  }

  // mult scalar
  d = DataVectorSP(d_rand);
  d.mult(scalar);
  for (int i = 0; i < N; ++i) {
    BOOST_CHECK_EQUAL(d[i], d_rand[i] * scalar);
  }

  // sum
  d = DataVectorSP(d_rand);
  BOOST_CHECK_CLOSE(d.sum(), sum, tol);

  // square
  d = DataVectorSP(d_rand);
  d.sqr();

  for (int i = 0; i < N; ++i) {
    BOOST_CHECK_EQUAL(d[i], d_rand[i] * d_rand[i]);
  }

  // sqrt
  d = DataVectorSP(d_rand);
  d.sqr();
  d.sqrt();
  for (int i = 0; i < N; ++i) {
    BOOST_CHECK_CLOSE(d[i], d_rand[i], tol);
  }

  // abs
  d = DataVectorSP(d_rand);
  d.abs();

  for (int i = 0; i < N; ++i) {
    BOOST_CHECK_EQUAL(d[i], std::abs(d_rand[i]));
  }

  // componentwise mult
  d = DataVectorSP(d_rand);

  d.componentwise_mult(d2);

  for (int i = 0; i < N; ++i) {
    BOOST_CHECK_EQUAL(d[i], d_rand[i] * static_cast<float>(i));
  }

  // componentwise div
  d = DataVectorSP(d_rand);

  for (int i = 0; i < N; ++i) {
    d2[i] = static_cast<float>(i) + 1.0f;
  }

  d.componentwise_div(d2);

  for (int i = 0; i < N; ++i) {
    BOOST_CHECK_CLOSE(d_rand[i] / (static_cast<float>(i) + 1.0f), d[i], tol);
  }
}

BOOST_AUTO_TEST_CASE(testDotProduct) {
  DataVectorSP d(3);
  float x = 0.0f;

  for (unsigned int i = 0; i < d.getSize(); ++i) {
    d[i] = static_cast<float>(i + 1);
    x += d[i] * d[i];
  }

  BOOST_CHECK_EQUAL(d.dotProduct(d), x);
}

BOOST_AUTO_TEST_SUITE_END()
