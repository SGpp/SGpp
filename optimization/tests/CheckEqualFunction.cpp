// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <sgpp/base/tools/RandomNumberGenerator.hpp>

#include <boost/test/unit_test.hpp>
#include <vector>

#include "CheckEqualFunction.hpp"

void checkEqualFunction(sgpp::base::ScalarFunction& f, sgpp::base::ScalarFunction& g) {
  sgpp::base::RandomNumberGenerator::getInstance().setSeed(42);

  const size_t d = f.getNumberOfParameters();
  const size_t N = 100;
  sgpp::base::DataVector x(d);

  BOOST_CHECK_EQUAL(g.getNumberOfParameters(), d);

  for (size_t i = 0; i < N; i++) {
    sgpp::base::RandomNumberGenerator::getInstance().getUniformRV(x);
    BOOST_CHECK_EQUAL(f.eval(x), g.eval(x));
  }
}

void checkEqualFunction(sgpp::base::ScalarFunctionGradient& f,
                        sgpp::base::ScalarFunctionGradient& g) {
  sgpp::base::RandomNumberGenerator::getInstance().setSeed(42);

  const size_t d = f.getNumberOfParameters();
  const size_t N = 100;
  sgpp::base::DataVector x(d);
  sgpp::base::DataVector gradient1(d), gradient2(d);

  BOOST_CHECK_EQUAL(g.getNumberOfParameters(), d);

  for (size_t i = 0; i < N; i++) {
    sgpp::base::RandomNumberGenerator::getInstance().getUniformRV(x);
    BOOST_CHECK_EQUAL(f.eval(x, gradient1), g.eval(x, gradient2));
    BOOST_CHECK_EQUAL(gradient2.getSize(), d);

    for (size_t t = 0; t < d; t++) {
      BOOST_CHECK_EQUAL(gradient1[t], gradient2[t]);
    }
  }
}

void checkEqualFunction(sgpp::base::ScalarFunctionHessian& f,
                        sgpp::base::ScalarFunctionHessian& g) {
  sgpp::base::RandomNumberGenerator::getInstance().setSeed(42);

  const size_t d = f.getNumberOfParameters();
  const size_t N = 100;
  sgpp::base::DataVector x(d);
  sgpp::base::DataVector gradient1(d), gradient2(d);
  sgpp::base::DataMatrix hessian1(d, d), hessian2(d, d);

  BOOST_CHECK_EQUAL(g.getNumberOfParameters(), d);

  for (size_t i = 0; i < N; i++) {
    sgpp::base::RandomNumberGenerator::getInstance().getUniformRV(x);
    BOOST_CHECK_EQUAL(f.eval(x, gradient1, hessian1), g.eval(x, gradient2, hessian2));
    BOOST_CHECK_EQUAL(gradient2.getSize(), d);
    BOOST_CHECK_EQUAL(hessian2.getNrows(), d);
    BOOST_CHECK_EQUAL(hessian2.getNcols(), d);

    for (size_t t1 = 0; t1 < d; t1++) {
      for (size_t t2 = 0; t2 < d; t2++) {
        BOOST_CHECK_EQUAL(hessian1(t1, t2), hessian2(t1, t2));
      }
    }
  }
}

void checkEqualFunction(sgpp::base::VectorFunction& f, sgpp::base::VectorFunction& g) {
  sgpp::base::RandomNumberGenerator::getInstance().setSeed(42);

  const size_t d = f.getNumberOfParameters();
  const size_t m = f.getNumberOfComponents();
  const size_t N = 100;
  sgpp::base::DataVector x(d);
  sgpp::base::DataVector value1(m), value2(m);

  BOOST_CHECK_EQUAL(g.getNumberOfParameters(), d);
  BOOST_CHECK_EQUAL(g.getNumberOfComponents(), m);

  for (size_t i = 0; i < N; i++) {
    sgpp::base::RandomNumberGenerator::getInstance().getUniformRV(x);
    f.eval(x, value1);
    g.eval(x, value2);

    BOOST_CHECK_EQUAL(value2.getSize(), m);

    for (size_t j = 0; j < m; j++) {
      BOOST_CHECK_EQUAL(value1[j], value2[j]);
    }
  }
}

void checkEqualFunction(sgpp::base::VectorFunctionGradient& f,
                        sgpp::base::VectorFunctionGradient& g) {
  sgpp::base::RandomNumberGenerator::getInstance().setSeed(42);

  const size_t d = f.getNumberOfParameters();
  const size_t m = f.getNumberOfComponents();
  const size_t N = 100;
  sgpp::base::DataVector x(d);
  sgpp::base::DataVector value1(m), value2(m);
  sgpp::base::DataMatrix gradient1(m, d), gradient2(m, d);

  BOOST_CHECK_EQUAL(g.getNumberOfParameters(), d);
  BOOST_CHECK_EQUAL(g.getNumberOfComponents(), m);

  for (size_t i = 0; i < N; i++) {
    sgpp::base::RandomNumberGenerator::getInstance().getUniformRV(x);
    f.eval(x, value1, gradient1);
    g.eval(x, value2, gradient2);

    BOOST_CHECK_EQUAL(value2.getSize(), m);
    BOOST_CHECK_EQUAL(gradient2.getNrows(), m);
    BOOST_CHECK_EQUAL(gradient2.getNcols(), d);

    for (size_t j = 0; j < m; j++) {
      BOOST_CHECK_EQUAL(value1[j], value2[j]);

      for (size_t t = 0; t < d; t++) {
        BOOST_CHECK_EQUAL(gradient1(j, t), gradient2(j, t));
      }
    }
  }
}

void checkEqualFunction(sgpp::base::VectorFunctionHessian& f,
                        sgpp::base::VectorFunctionHessian& g) {
  sgpp::base::RandomNumberGenerator::getInstance().setSeed(42);

  const size_t d = f.getNumberOfParameters();
  const size_t m = f.getNumberOfComponents();
  const size_t N = 100;
  sgpp::base::DataVector x(d);
  sgpp::base::DataVector value1(m), value2(m);
  sgpp::base::DataMatrix gradient1(m, d), gradient2(m, d);
  std::vector<sgpp::base::DataMatrix> hessian1(m, sgpp::base::DataMatrix(d, d));
  std::vector<sgpp::base::DataMatrix> hessian2(m, sgpp::base::DataMatrix(d, d));

  BOOST_CHECK_EQUAL(g.getNumberOfParameters(), d);
  BOOST_CHECK_EQUAL(g.getNumberOfComponents(), m);

  for (size_t i = 0; i < N; i++) {
    sgpp::base::RandomNumberGenerator::getInstance().getUniformRV(x);
    f.eval(x, value1, gradient1, hessian1);
    g.eval(x, value2, gradient2, hessian2);

    BOOST_CHECK_EQUAL(value2.getSize(), m);
    BOOST_CHECK_EQUAL(gradient2.getNrows(), m);
    BOOST_CHECK_EQUAL(gradient2.getNcols(), d);
    BOOST_CHECK_EQUAL(hessian2.size(), m);

    for (size_t j = 0; j < m; j++) {
      BOOST_CHECK_EQUAL(value1[j], value2[j]);
      BOOST_CHECK_EQUAL(hessian2[j].getNrows(), d);
      BOOST_CHECK_EQUAL(hessian2[j].getNcols(), d);

      for (size_t t = 0; t < d; t++) {
        BOOST_CHECK_EQUAL(gradient1(j, t), gradient2(j, t));

        for (size_t t2 = 0; t2 < d; t2++) {
          BOOST_CHECK_EQUAL(hessian1[j](t, t2), hessian2[j](t, t2));
        }
      }
    }
  }
}
