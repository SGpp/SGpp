#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
using namespace SGPP::base;

struct FixtureDataMatrix {
	FixtureDataMatrix() :
    nrows(5), ncols(3), N(nrows * ncols), d_rand(nrows, ncols), min(0), max(0), sum(
      0) {
    //l_rand_total = new double[nrows * ncols];
    l_rand = new double* [nrows];

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
  ~FixtureDataMatrix() {

    for (int i = 0; i < nrows; ++i) {
      delete [] l_rand[i];
    }

    delete [] l_rand;
    BOOST_TEST_MESSAGE("teardown fixture");
  }
  int nrows, ncols, N;
  double** l_rand;
  DataMatrix d_rand;
  double min, max, sum;
};

BOOST_FIXTURE_TEST_SUITE(testDataMatrix, FixtureDataMatrix)

BOOST_AUTO_TEST_CASE(testConstructor) {
  DataMatrix d(42,17);
  BOOST_CHECK_EQUAL(d.getSize(), 42*17);
  BOOST_CHECK_EQUAL(d.getNrows(), 42);
  BOOST_CHECK_EQUAL(d.getNcols(), 17);
}

BOOST_AUTO_TEST_CASE(testSetUp) {
  for (int i = 0; i < nrows; ++i) {
    for (int j = 0; j < ncols; ++j) {
	  BOOST_CHECK_EQUAL(d_rand.get(i,j), l_rand[i][j]);
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

	for (int i = 0; i < nrows; ++i) {
	  for (int j = 0; j < ncols; ++j) {
	    d2.set(i,j, rand()+static_cast<double>(i*j));
	  }
    }

	//abs
	d.abs();
	for (int i = 0; i < nrows; ++i) {
	  for (int j = 0; j < ncols; ++j) {
		  BOOST_CHECK_EQUAL(fabs(d_rand.get(i,j)), d.get(i,j));
	  }
    }

	//add
	d = DataMatrix(d_rand);
	d.add(d2);
	for (int i = 0; i < nrows; ++i) {
	  for (int j = 0; j < ncols; ++j) {
		  BOOST_CHECK_EQUAL(d_rand.get(i,j)+d2.get(i,j), d.get(i,j));
	  }
    }

	//addReduce
	d = DataMatrix(d_rand);
	DataVector reduction(nrows);
	d.addReduce(reduction);
	double reduce_sum = 0.0;
	for (int i = 0; i < nrows; ++i) {
	  reduce_sum = 0.0;
	  for (int j = 0; j < ncols; ++j) {
		  reduce_sum += d_rand.get(i,j);
	  }
	  BOOST_CHECK_CLOSE(reduce_sum, reduction[i], tol);
	}

	//componentwise_div
	d = DataMatrix(d_rand);
	d.componentwise_div(d2);
	for (int i = 0; i < nrows; ++i) {
	  for (int j = 0; j < ncols; ++j) {
		  BOOST_CHECK_CLOSE( d_rand.get(i,j)*(1.0/d2.get(i,j)), d.get(i,j), tol );
	  }
	}

	//componentwise_mult
	d = DataMatrix(d_rand);
	d.componentwise_mult(d2);
	for (int i = 0; i < nrows; ++i){
		for (int j = 0; j < ncols; ++j){
			BOOST_CHECK_EQUAL(d_rand.get(i,j)*d2.get(i,j), d.get(i,j));
		}
	}

	//max column
	d = DataMatrix(d_rand);
	double max = 0.0;
	for (int j = 0; j < ncols; ++j){
		max = 0.0;
		for (int i = 0; i < nrows; ++i){
			max = max < d_rand.get(i,j) ? d_rand.get(i,j) : max;
		}
		BOOST_CHECK_EQUAL(max, d.max(j));
	}

}

BOOST_AUTO_TEST_SUITE_END()

