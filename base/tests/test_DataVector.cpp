#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp/base/datatypes/DataVector.hpp>
using namespace SGPP::base;

struct Fixture {
  Fixture() :
    nrows(5), ncols(4), N(nrows * ncols), d_rand(N), min(0), max(0), sum(
      0) {
    l_rand_total = new SGPP::float_t[nrows * ncols];
    l_rand = new SGPP::float_t* [nrows];

    for (int i = 0; i < nrows; ++i) {
      l_rand[i] = new SGPP::float_t[ncols];
    }

    for (int i = 0; i < nrows; ++i) {
      for (int j = 0; j < ncols; ++j) {
        l_rand[i][j] = i * j + i * 0.5 + 2.34 * j;
        l_rand_total[i * ncols + j] = l_rand[i][j];
        min = min > l_rand[i][j] ? l_rand[i][j] : min;
        max = max < l_rand[i][j] ? l_rand[i][j] : max;
        sum += l_rand[i][j];
      }
    }

    for (int i = 0; i < N; ++i) {
      d_rand[i] = l_rand_total[i];
    }

    BOOST_TEST_MESSAGE("setup fixture");
  }
  ~Fixture() {
    delete l_rand_total;

    for (int i = 0; i < nrows; ++i) {
      delete l_rand[i];
    }

    delete l_rand;
    BOOST_TEST_MESSAGE("teardown fixture");
  }
  int nrows, ncols, N;
  SGPP::float_t** l_rand;
  SGPP::float_t* l_rand_total;
  DataVector d_rand;
  SGPP::float_t min, max, sum;
};

BOOST_FIXTURE_TEST_SUITE(testDataVector, Fixture)

BOOST_AUTO_TEST_CASE(testConstructor) {
  DataVector d = DataVector(2);
  BOOST_CHECK_EQUAL(d.getSize(), 2U);
}

BOOST_AUTO_TEST_CASE(testSetUp) {
  for (int i = 0; i < N; ++i) {
    BOOST_CHECK_EQUAL(d_rand[i], l_rand_total[i]);
  }
}
BOOST_AUTO_TEST_CASE(testMinMax) {
  BOOST_CHECK_EQUAL(d_rand.min(), min);
  BOOST_CHECK_EQUAL(d_rand.max(), max);
}

BOOST_AUTO_TEST_CASE(testOps) {
  //sum
  BOOST_CHECK_EQUAL(d_rand.sum(), sum);
  //square
  DataVector d = d_rand;
  d.sqr();

  for (int i = 0; i < N; ++i) {
    BOOST_CHECK_EQUAL(d_rand[i] * d_rand[i], d[i]);
  }

  //abs
  d = DataVector(d_rand);
  d.abs();

  for (int i = 0; i < N; ++i) {
    BOOST_CHECK_EQUAL(fabs(d_rand[i]), d[i]);
  }

  //componentwise mult
  d = DataVector(d_rand);
  DataVector d2 = DataVector(N);

  for (int i = 0; i < N; ++i) {
    d2[i] = static_cast<SGPP::float_t>(i);
  }

  d.componentwise_mult(d2);

  for (int i = 0; i < N; ++i) {
    BOOST_CHECK_EQUAL(d_rand[i] * static_cast<SGPP::float_t>(i), d[i]);
  }

  // componentwise div
  d = DataVector(d_rand);

  for (int i = 0; i < N; ++i) {
    d2[i] = i + 1.0;
  }

  d.componentwise_div(d2);

  for (int i = 0; i < N; ++i) {
#if USE_DOUBLE_PRECISION == 1
    BOOST_CHECK_CLOSE( d_rand[i] / (i + 1.0), d[i], float_t(1e-12) );
#else
    BOOST_CHECK_CLOSE( d_rand[i] / (i + 1.0), d[i], float_t(1e-5) );
#endif
  }
}

BOOST_AUTO_TEST_CASE(testDotProduct) {
  DataVector d = DataVector(3);
  SGPP::float_t x = 0;

  for (unsigned int i = 0; i < d.getSize(); ++i) {
    d[i] = static_cast<SGPP::float_t>(i + 1);
    x += d[i] * d[i];
  }

  BOOST_CHECK_EQUAL(d.dotProduct(d), x);
}

BOOST_AUTO_TEST_SUITE_END()

