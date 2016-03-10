// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp/base/datatypes/DataMatrixSP.hpp>
#include <sgpp/base/datatypes/DataVectorSP.hpp>

#include <algorithm>

using sgpp::base::DataMatrixSP;
using sgpp::base::DataVectorSP;

struct FixtureDataMatrixSP {
  FixtureDataMatrixSP()
      : nrows(5), ncols(3), N(nrows * ncols), d_rand(nrows, ncols), min(0), max(0), sum(0) {
    l_rand = new float*[nrows];

    for (int i = 0; i < nrows; ++i) {
      l_rand[i] = new float[ncols];
    }

    float value;
    for (int i = 0; i < nrows; ++i) {
      for (int j = 0; j < ncols; ++j) {
        value = static_cast<float>(i) * static_cast<float>(j) + static_cast<float>(i) * 0.5f +
                2.34f * static_cast<float>(j);
        l_rand[i][j] = value;
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
  ~FixtureDataMatrixSP() {
    for (int i = 0; i < nrows; ++i) {
      delete[] l_rand[i];
    }

    delete[] l_rand;
    BOOST_TEST_MESSAGE("teardown fixture");
  }
  int nrows, ncols, N;
  float** l_rand;
  DataMatrixSP d_rand;
  float min, max, sum;
};

BOOST_FIXTURE_TEST_SUITE(testDataMatrixSP, FixtureDataMatrixSP)

BOOST_AUTO_TEST_CASE(testConstructor) {
  DataMatrixSP d(42, 17);
  BOOST_CHECK_EQUAL(d.getSize(), 42 * 17);
  BOOST_CHECK_EQUAL(d.getNrows(), 42);
  BOOST_CHECK_EQUAL(d.getNcols(), 17);
}

BOOST_AUTO_TEST_CASE(testSetUp) {
  for (int i = 0; i < nrows; ++i) {
    for (int j = 0; j < ncols; ++j) {
      BOOST_CHECK_EQUAL(d_rand.get(i, j), l_rand[i][j]);
    }
  }
}

BOOST_AUTO_TEST_CASE(testMinMax) {
  BOOST_CHECK_EQUAL(d_rand.min(), min);
  BOOST_CHECK_EQUAL(d_rand.max(), max);
}

BOOST_AUTO_TEST_CASE(testOps) {
  float tol = 0.00002f;

  DataMatrixSP d = d_rand;
  DataMatrixSP d2(nrows, ncols);

  for (int i = 0; i < nrows; ++i) {
    for (int j = 0; j < ncols; ++j) {
      d2.set(i, j, 1.0f + 2.123f * (static_cast<float>(rand()) / static_cast<float>(RAND_MAX)) +
                       static_cast<float>(i * j));
    }
  }

  // abs
  d.abs();
  for (int i = 0; i < nrows; ++i) {
    for (int j = 0; j < ncols; ++j) {
      BOOST_CHECK_EQUAL(d.get(i, j), fabs(d_rand.get(i, j)));
    }
  }

  // add
  d = DataMatrixSP(d_rand);
  d.add(d2);
  for (int i = 0; i < nrows; ++i) {
    for (int j = 0; j < ncols; ++j) {
      BOOST_CHECK_EQUAL(d.get(i, j), d_rand.get(i, j) + d2.get(i, j));
    }
  }

  // addReduce
  d = DataMatrixSP(d_rand);
  DataVectorSP reduction(nrows);
  d.addReduce(reduction);
  float reduce_sum = 0.0;
  for (int i = 0; i < nrows; ++i) {
    reduce_sum = 0.0;
    for (int j = 0; j < ncols; ++j) {
      reduce_sum += d_rand.get(i, j);
    }
    BOOST_CHECK_CLOSE(reduction[i], reduce_sum, tol);
  }

  // componentwise_div
  d = DataMatrixSP(d_rand);
  d.componentwise_div(d2);
  for (int i = 0; i < nrows; ++i) {
    for (int j = 0; j < ncols; ++j) {
      BOOST_CHECK_CLOSE(d.get(i, j), d_rand.get(i, j) * (1.0 / d2.get(i, j)), tol);
    }
  }

  // componentwise_mult
  d = DataMatrixSP(d_rand);
  d.componentwise_mult(d2);
  for (int i = 0; i < nrows; ++i) {
    for (int j = 0; j < ncols; ++j) {
      BOOST_CHECK_EQUAL(d.get(i, j), d_rand.get(i, j) * d2.get(i, j));
    }
  }

  // max column
  d = DataMatrixSP(d_rand);
  float colMax;
  for (int j = 0; j < ncols; ++j) {
    colMax = d.get(0, j);
    for (int i = 0; i < nrows; ++i) {
      colMax = colMax < d_rand.get(i, j) ? d_rand.get(i, j) : colMax;
    }
    BOOST_CHECK_EQUAL(d.max(j), colMax);
  }

  // total max
  d = DataMatrixSP(d_rand);
  float d_max = d.max();
  BOOST_CHECK_EQUAL(d_max, max);

  // min column
  d = DataMatrixSP(d_rand);
  float colMin;
  for (int j = 0; j < ncols; ++j) {
    colMin = d.get(0, j);
    for (int i = 0; i < nrows; ++i) {
      colMin = colMin > d_rand.get(i, j) ? d_rand.get(i, j) : colMin;
    }
    BOOST_CHECK_EQUAL(d.min(j), colMin);
  }

  // total min
  d = DataMatrixSP(d_rand);
  float d_min = d.min();
  BOOST_CHECK_EQUAL(d_min, min);

  // min max column
  float minActual = 0;
  float maxActual = 0;
  for (int j = 0; j < ncols; ++j) {
    colMax = d.get(0, j);
    colMin = d.get(0, j);
    d.minmax(j, &minActual, &maxActual);
    for (int i = 0; i < nrows; ++i) {
      colMin = colMin > d_rand.get(i, j) ? d_rand.get(i, j) : colMin;
      colMax = colMax < d_rand.get(i, j) ? d_rand.get(i, j) : colMax;
    }
    BOOST_CHECK_EQUAL(minActual, colMin);
    BOOST_CHECK_EQUAL(maxActual, colMax);
  }

  // min max all
  d.minmax(&minActual, &maxActual);
  BOOST_CHECK_EQUAL(minActual, min);
  BOOST_CHECK_EQUAL(maxActual, max);

  // mult scalar
  d = DataMatrixSP(d_rand);
  float scalar = 1.3124f;
  d.mult(scalar);
  for (int i = 0; i < nrows; ++i) {
    for (int j = 0; j < ncols; ++j) {
      BOOST_CHECK_CLOSE(d.get(i, j), d_rand.get(i, j) * scalar, tol);
    }
  }

  // normalize dimension of column j
  d = DataMatrixSP(d_rand);
  float border = 0.0f;
  float delta;
  for (int j = 0; j < ncols; ++j) {
    d.normalizeDimension(j);
    delta = (d_rand.max(j) - d_rand.min(j)) / (1 - 2 * border);
    for (int i = 0; i < nrows; ++i) {
      BOOST_CHECK_CLOSE(d.get(i, j), (d_rand.get(i, j) - d_rand.min(j)) / delta + border, tol);
    }
  }

  // normalize dimension with border
  d = DataMatrixSP(d_rand);
  border = 3.1415f;
  for (int j = 0; j < ncols; ++j) {
    d.normalizeDimension(j, border);
    delta = (d_rand.max(j) - d_rand.min(j)) / (1 - 2 * border);
    for (int i = 0; i < nrows; ++i) {
      BOOST_CHECK_CLOSE(d.get(i, j), (d_rand.get(i, j) - d_rand.min(j)) / delta + border, tol);
    }
  }

  // resize rows
  d = DataMatrixSP(d_rand);
  const int newRows = nrows - 1 < 1 ? 1 : nrows - 1;
  d.resize(newRows);
  BOOST_CHECK_EQUAL(d.getNrows(), newRows);

  // resizeZero rows
  d.resize(nrows, ncols);
  d = DataMatrixSP(d_rand);
  d.resizeZero(newRows);
  BOOST_CHECK_EQUAL(d.getNrows(), newRows);

  // resize rows and cols
  d.resize(nrows, ncols);
  d = DataMatrixSP(d_rand);
  const int newCols = ncols - 1 < 1 ? 1 : ncols - 1;
  d.resize(newRows, newCols);
  BOOST_CHECK_EQUAL(d.getNrows(), newRows);
  BOOST_CHECK_EQUAL(d.getNcols(), newCols);

  // resizeZero rows and cols
  d.resize(nrows, ncols);
  d = DataMatrixSP(d_rand);
  d.resizeZero(newRows, newCols);
  BOOST_CHECK_EQUAL(d.getNrows(), newRows);
  BOOST_CHECK_EQUAL(d.getNcols(), newCols);
  d.resize(nrows, ncols);

  // setAll
  d = DataMatrixSP(d_rand);
  float setValue = 3.12f;
  d.setAll(setValue);
  for (int i = 0; i < nrows; ++i) {
    for (int j = 0; j < ncols; ++j) {
      BOOST_CHECK_EQUAL(d.get(i, j), setValue);
    }
  }

  // setColumn
  d = DataMatrixSP(d_rand);
  DataVectorSP setVectorRow(nrows);
  setVectorRow.setAll(setValue + 1);
  for (int j = 0; j < ncols; ++j) {
    d.setColumn(j, setVectorRow);
    for (int i = 0; i < nrows; ++i) {
      BOOST_CHECK_EQUAL(d.get(i, j), setValue + 1);
    }
  }

  // setRow
  d = DataMatrixSP(d_rand);
  DataVectorSP setVectorCol(ncols);
  setVectorCol.setAll(setValue + 1);
  for (int i = 0; i < nrows; ++i) {
    d.setRow(i, setVectorCol);
    for (int j = 0; j < ncols; ++j) {
      BOOST_CHECK_EQUAL(d.get(i, j), setValue + 1);
    }
  }

  // sqr
  d = DataMatrixSP(d_rand);
  d.sqr();
  for (int i = 0; i < nrows; ++i) {
    for (int j = 0; j < ncols; ++j) {
      BOOST_CHECK_CLOSE(d.get(i, j), d_rand.get(i, j) * d_rand.get(i, j), tol);
    }
  }

  // sqrt
  d = DataMatrixSP(d_rand);
  d.sqrt();
  d.sqr();
  float sqrtTol = 0.0001f;
  for (int i = 0; i < nrows; ++i) {
    for (int j = 0; j < ncols; ++j) {
      BOOST_CHECK_CLOSE(d.get(i, j), d_rand.get(i, j), sqrtTol);
    }
  }

  // sub
  d = DataMatrixSP(d_rand);
  d.sub(d2);
  for (int i = 0; i < nrows; ++i) {
    for (int j = 0; j < ncols; ++j) {
      BOOST_CHECK_EQUAL(d.get(i, j), d_rand.get(i, j) - d2.get(i, j));
    }
  }

  // sum
  d = DataMatrixSP(d_rand);
  float actualSum = d.sum();
  float expectedSum = 0.0f;
  for (int i = 0; i < nrows; ++i) {
    for (int j = 0; j < ncols; ++j) {
      expectedSum += d_rand.get(i, j);
    }
  }
  BOOST_CHECK_CLOSE(actualSum, expectedSum, tol);

  // transpose
  d = DataMatrixSP(d_rand);
  d.transpose();
  for (int i = 0; i < nrows; ++i) {
    for (int j = 0; j < ncols; ++j) {
      BOOST_CHECK_EQUAL(d.get(j, i), d_rand(i, j));
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
