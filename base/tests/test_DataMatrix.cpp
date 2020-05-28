// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <algorithm>
#include <cmath>

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;

struct FixtureDataMatrix {
  FixtureDataMatrix()
      : nrows(5), ncols(3), N(nrows * ncols), d_rand(nrows, ncols), min(0), max(0), sum(0) {
    l_rand = new double*[nrows];

    for (int i = 0; i < nrows; i++) {
      l_rand[i] = new double[ncols];
    }

    for (int i = 0; i < nrows; i++) {
      for (int j = 0; j < ncols; j++) {
        l_rand[i][j] = static_cast<double>(i) * static_cast<double>(j) +
                static_cast<double>(i) * 0.5 + 2.34 * static_cast<double>(j);
        min = min > l_rand[i][j] ? l_rand[i][j] : min;
        max = max < l_rand[i][j] ? l_rand[i][j] : max;
        sum += l_rand[i][j];
      }
    }

    for (int i = 0; i < nrows; i++) {
      for (int j = 0; j < ncols; j++) {
        d_rand.set(i, j, l_rand[i][j]);
      }
    }

    BOOST_TEST_MESSAGE("setup fixture");
  }
  ~FixtureDataMatrix() {
    for (int i = 0; i < nrows; i++) {
      delete[] l_rand[i];
    }

    delete[] l_rand;
    BOOST_TEST_MESSAGE("teardown fixture");
  }
  int nrows, ncols, N;
  double** l_rand;
  DataMatrix d_rand;
  double min, max, sum;
};

BOOST_FIXTURE_TEST_SUITE(testDataMatrix, FixtureDataMatrix)

BOOST_AUTO_TEST_CASE(testConstructor) {
  DataMatrix d(42, 17);
  BOOST_CHECK_EQUAL(d.getSize(), 42 * 17);
  BOOST_CHECK_EQUAL(d.getNrows(), 42);
  BOOST_CHECK_EQUAL(d.getNcols(), 17);

  d = DataMatrix({1.0, 0.5, -2.0, 2.5, 1.0, -0.5}, 3);
  BOOST_CHECK_EQUAL(d.getSize(), 6);
  BOOST_CHECK_EQUAL(d.getNrows(), 3);
  BOOST_CHECK_EQUAL(d.getNcols(), 2);
  BOOST_CHECK_EQUAL(d(0, 0), 1.0);
  BOOST_CHECK_EQUAL(d(0, 1), 0.5);
  BOOST_CHECK_EQUAL(d(1, 0), -2.0);
  BOOST_CHECK_EQUAL(d(1, 1), 2.5);
  BOOST_CHECK_EQUAL(d(2, 0), 1.0);
  BOOST_CHECK_EQUAL(d(2, 1), -0.5);
}

BOOST_AUTO_TEST_CASE(testSetUp) {
  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < ncols; j++) {
      BOOST_CHECK_EQUAL(d_rand.get(i, j), l_rand[i][j]);
    }
  }
}

BOOST_AUTO_TEST_CASE(testMinMax) {
  BOOST_CHECK_EQUAL(d_rand.min(), min);
  BOOST_CHECK_EQUAL(d_rand.max(), max);
}

BOOST_AUTO_TEST_CASE(testOps) {
  double tol = 1e-12;

  DataMatrix d = d_rand;
  DataMatrix d2(nrows, ncols);

  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < ncols; j++) {
      d2.set(i, j, 1.0 + 2.123 * (static_cast<double>(rand()) / static_cast<double>(RAND_MAX)) +
                       static_cast<double>(i * j));
    }
  }

  // abs
  d.abs();
  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < ncols; j++) {
      BOOST_CHECK_EQUAL(d.get(i, j), std::abs(d_rand.get(i, j)));
    }
  }

  // add
  d = DataMatrix(d_rand);
  d.add(d2);
  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < ncols; j++) {
      BOOST_CHECK_EQUAL(d.get(i, j), d_rand.get(i, j) + d2.get(i, j));
    }
  }

  // addReduce
  d = DataMatrix(d_rand);
  DataVector reduction(nrows);
  d.addReduce(reduction);
  double reduce_sum = 0.0;
  for (int i = 0; i < nrows; i++) {
    reduce_sum = 0.0;
    for (int j = 0; j < ncols; j++) {
      reduce_sum += d_rand.get(i, j);
    }
    BOOST_CHECK_CLOSE(reduction[i], reduce_sum, tol);
  }

  // componentwise_div
  d = DataMatrix(d_rand);
  d.componentwise_div(d2);
  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < ncols; j++) {
      BOOST_CHECK_CLOSE(d.get(i, j), d_rand.get(i, j) * (1.0 / d2.get(i, j)), tol);
    }
  }

  // componentwise_mult
  d = DataMatrix(d_rand);
  d.componentwise_mult(d2);
  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < ncols; j++) {
      BOOST_CHECK_EQUAL(d.get(i, j), d_rand.get(i, j) * d2.get(i, j));
    }
  }

  // max column
  d = DataMatrix(d_rand);
  double colMax;
  for (int j = 0; j < ncols; j++) {
    colMax = d.get(0, j);
    for (int i = 0; i < nrows; i++) {
      colMax = colMax < d_rand.get(i, j) ? d_rand.get(i, j) : colMax;
    }
    BOOST_CHECK_EQUAL(d.max(j), colMax);
  }

  // total max
  d = DataMatrix(d_rand);
  double d_max = d.max();
  BOOST_CHECK_EQUAL(d_max, max);

  // min column
  d = DataMatrix(d_rand);
  double colMin;
  for (int j = 0; j < ncols; j++) {
    colMin = d.get(0, j);
    for (int i = 0; i < nrows; i++) {
      colMin = colMin > d_rand.get(i, j) ? d_rand.get(i, j) : colMin;
    }
    BOOST_CHECK_EQUAL(d.min(j), colMin);
  }

  // total min
  d = DataMatrix(d_rand);
  double d_min = d.min();
  BOOST_CHECK_EQUAL(d_min, min);

  // min max column
  double minActual = 0.0;
  double maxActual = 0.0;
  for (int j = 0; j < ncols; j++) {
    colMax = d.get(0, j);
    colMin = d.get(0, j);
    d.minmax(j, &minActual, &maxActual);
    for (int i = 0; i < nrows; i++) {
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
  d = DataMatrix(d_rand);
  double scalar = 1.3124;
  d.mult(scalar);
  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < ncols; j++) {
      BOOST_CHECK_CLOSE(d.get(i, j), d_rand.get(i, j) * scalar, tol);
    }
  }

  // normalize dimension of column j
  d = DataMatrix(d_rand);
  double border = 0.0;
  double delta;
  for (int j = 0; j < ncols; j++) {
    d.normalizeDimension(j);
    delta = (d_rand.max(j) - d_rand.min(j)) / (1.0 - 2.0 * border);
    for (int i = 0; i < nrows; i++) {
      BOOST_CHECK_CLOSE(d.get(i, j), (d_rand.get(i, j) - d_rand.min(j)) / delta + border, tol);
    }
  }

  // normalize dimension with border
  d = DataMatrix(d_rand);
  border = 3.1415;
  for (int j = 0; j < ncols; j++) {
    d.normalizeDimension(j, border);
    delta = (d_rand.max(j) - d_rand.min(j)) / (1.0 - 2.0 * border);
    for (int i = 0; i < nrows; i++) {
      BOOST_CHECK_CLOSE(d.get(i, j), (d_rand.get(i, j) - d_rand.min(j)) / delta + border, tol);
    }
  }

  // resize rows
  d = DataMatrix(d_rand);
  const int newRows = nrows - 1 < 1 ? 1 : nrows - 1;
  d.resize(newRows);
  BOOST_CHECK_EQUAL(d.getNrows(), newRows);

  // resizeZero rows
  d.resize(nrows, ncols);
  d = DataMatrix(d_rand);
  d.resizeZero(newRows);
  BOOST_CHECK_EQUAL(d.getNrows(), newRows);

  // resize rows and cols
  d.resize(nrows, ncols);
  d = DataMatrix(d_rand);
  const int newCols = ncols - 1 < 1 ? 1 : ncols - 1;
  d.resize(newRows, newCols);
  BOOST_CHECK_EQUAL(d.getNrows(), newRows);
  BOOST_CHECK_EQUAL(d.getNcols(), newCols);

  // resizeZero rows and cols
  d.resize(nrows, ncols);
  d = DataMatrix(d_rand);
  d.resizeZero(newRows, newCols);
  BOOST_CHECK_EQUAL(d.getNrows(), newRows);
  BOOST_CHECK_EQUAL(d.getNcols(), newCols);
  d.resize(nrows, ncols);

  // setAll
  d = DataMatrix(d_rand);
  double setValue = 3.12;
  d.setAll(setValue);
  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < ncols; j++) {
      BOOST_CHECK_EQUAL(d.get(i, j), setValue);
    }
  }

  // setColumn
  d = DataMatrix(d_rand);
  DataVector setVectorRow(nrows);
  setVectorRow.setAll(setValue + 1.0);
  for (int j = 0; j < ncols; j++) {
    d.setColumn(j, setVectorRow);
    for (int i = 0; i < nrows; i++) {
      BOOST_CHECK_EQUAL(d.get(i, j), setValue + 1.0);
    }
  }

  // setRow
  d = DataMatrix(d_rand);
  DataVector setVectorCol(ncols);
  setVectorCol.setAll(setValue + 1.0);
  for (int i = 0; i < nrows; i++) {
    d.setRow(i, setVectorCol);
    for (int j = 0; j < ncols; j++) {
      BOOST_CHECK_EQUAL(d.get(i, j), setValue + 1.0);
    }
  }

  // sqr
  d = DataMatrix(d_rand);
  d.sqr();
  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < ncols; j++) {
      BOOST_CHECK_CLOSE(d.get(i, j), d_rand.get(i, j) * d_rand.get(i, j), tol);
    }
  }

  // sqrt
  d = DataMatrix(d_rand);
  d.sqrt();
  d.sqr();
  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < ncols; j++) {
      BOOST_CHECK_CLOSE(d.get(i, j), d_rand.get(i, j), tol);
    }
  }

  // sub
  d = DataMatrix(d_rand);
  d.sub(d2);
  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < ncols; j++) {
      BOOST_CHECK_EQUAL(d.get(i, j), d_rand.get(i, j) - d2.get(i, j));
    }
  }

  // sum
  d = DataMatrix(d_rand);
  double actualSum = d.sum();
  double expectedSum = 0.0;
  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < ncols; j++) {
      expectedSum += d_rand.get(i, j);
    }
  }
  BOOST_CHECK_CLOSE(actualSum, expectedSum, tol);

  // transpose
  d = DataMatrix(d_rand);
  d.transpose();
  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < ncols; j++) {
      BOOST_CHECK_EQUAL(d.get(j, i), d_rand(i, j));
    }
  }
}

BOOST_AUTO_TEST_CASE(quadraticTransposeTest) {
  for (size_t rows : {0, 1, 2, 3, 5, 7, 10, 20}) {
    for (size_t cols : {0, 1, 2, 3, 5, 7, 10, 20}) {
      DataMatrix m(rows, cols);
      {
        size_t k = 0;
        for (size_t i = 0; i < rows; i++) {
          for (size_t j = 0; j < cols; j++) {
            m(i, j) = static_cast<double>(k);
            k++;
          }
        }
      }
      DataMatrix t = m;
      t.transpose();
      for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
          BOOST_CHECK_MESSAGE(m(i, j) == t(j, i), "rows=" << rows << " cols=" << cols << " m(" << i
                                                          << "," << j << ")=" << m(i, j)
                                                          << "!=" << t(j, i));
        }
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(appendToColTest) {
  for (size_t rows = 0; rows < 20; rows++) {
    for (size_t cols = 0; cols < 20; cols++) {
      DataMatrix m(rows, cols);
      DataVector v(rows);
      {
        size_t k = 0;
        for (size_t i = 0; i < rows; i++) {
          for (size_t j = 0; j < cols; j++) {
            m(i, j) = static_cast<double>(k);
            k++;
          }
          v[i] = static_cast<double>(k);
          k++;
        }
      }
      m.appendCol(v);
      {
        size_t k = 0;
        for (size_t i = 0; i < rows; i++) {
          for (size_t j = 0; j < cols + 1; j++) {
            BOOST_CHECK_MESSAGE(m(i, j) == static_cast<double>(k),
                                "rows=" << rows << " cols=" << cols << " m(" << i
                                        << "," << j << ")=" << m(i, j) << "!=" << k);
            k++;
          }
        }
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(resizeToSubMatrixTests) {
  size_t rows = 20;
  size_t cols = 10;

  DataMatrix m(rows, cols);
  {
    size_t k = 0;
    for (size_t i = 0; i < rows; i++) {
      for (size_t j = 0; j < cols; j++) {
        m(i, j) = static_cast<double>(k);
        k++;
      }
    }
  }
  size_t x1 = 12;
  size_t x2 = 18;
  size_t y1 = 4;
  size_t y2 = 9;
  m.resizeToSubMatrix(x1, y1, x2, y2);
  BOOST_CHECK_EQUAL(m.getNrows(), x2 - x1 + 1);
  BOOST_CHECK_EQUAL(m.getNcols(), y2 - y1 + 1);
  {
    for (size_t i = 0; i < x2 - x1 + 1; i++) {
      for (size_t j = 0; j < y2 - y1 + 1; j++) {
        BOOST_CHECK_EQUAL(m(i, j), (i + x1) * cols + y1 + j);
      }
    }
  }

  DataMatrix m_narrow(3, 1);
  {
    m_narrow(0, 0) = 0.0;
    m_narrow(1, 0) = 1.0;
    m_narrow(2, 0) = 2.0;
  }
  m_narrow.resizeToSubMatrix(1, 0, 2, 0);
  {
    BOOST_CHECK_EQUAL(m_narrow(0, 0), 1.0);
    BOOST_CHECK_EQUAL(m_narrow(1, 0), 2.0);
  }
}

BOOST_AUTO_TEST_SUITE_END()
