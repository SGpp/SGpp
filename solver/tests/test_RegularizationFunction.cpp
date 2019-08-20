// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/solver/sle/fista/ElasticNetFunction.hpp>
#include <sgpp/solver/sle/fista/LassoFunction.hpp>
#include <sgpp/solver/sle/fista/RidgeFunction.hpp>
#include <sgpp/solver/sle/fista/ZeroFunction.hpp>

using sgpp::base::DataVector;
using sgpp::solver::ZeroFunction;
using sgpp::solver::ElasticNetFunction;
using sgpp::solver::LassoFunction;
using sgpp::solver::RidgeFunction;

BOOST_AUTO_TEST_SUITE(TestRegularizationFunction)

const double eps = 1e-5;
const double stepsize = 1.0;
const double lambda = 1.0;
const auto testVec = DataVector(5, -2.0);

BOOST_AUTO_TEST_CASE(testZeroFunction) {
  auto fun = ZeroFunction();
  const double eval = fun.eval(testVec);
  const auto prox = fun.prox(testVec, stepsize);
  BOOST_CHECK_CLOSE(eval, 0.0, eps);

  for (size_t i = 0; i < testVec.getSize(); ++i) {
    BOOST_CHECK_CLOSE(prox[i], (-2.0), eps);
  }
}

BOOST_AUTO_TEST_CASE(test_ElasticNetFunction) {
  auto enRidge = ElasticNetFunction(lambda, 0.0);
  auto enLasso = ElasticNetFunction(lambda, 1.0);
  auto enZero = ElasticNetFunction(0.0, 0.5);
  auto enMixed = ElasticNetFunction(lambda, 0.5);

  auto ridge = RidgeFunction(lambda);
  auto lasso = LassoFunction(lambda);
  auto zero = ZeroFunction();

  const double evalEnRidge = enRidge.eval(testVec);
  const double evalRidge = ridge.eval(testVec);
  const double evalEnLasso = enLasso.eval(testVec);
  const double evalEnMixed = enMixed.eval(testVec);
  const double evalLasso = lasso.eval(testVec);
  const double evalEnZero = enZero.eval(testVec);
  const double evalZero = zero.eval(testVec);
  const double evalMixed = 0.5 * (20) + 0.5 * (10);

  const auto proxEnRidge = enRidge.prox(testVec, stepsize);
  const auto proxRidge = ridge.prox(testVec, stepsize);
  const auto proxEnLasso = enLasso.prox(testVec, stepsize);
  const auto proxLasso = lasso.prox(testVec, stepsize);
  const auto proxEnZero = enZero.prox(testVec, stepsize);
  const auto proxZero = zero.prox(testVec, stepsize);
  const auto proxEnMixed = enMixed.prox(testVec, stepsize);
  const auto proxMixed = DataVector(5, 0.5 * (-1.5));
  // First ridge: 0.5 (with lambda = 0.5)
  // Then lasso: -1.5 (with lambda = 0.5)

  BOOST_CHECK_CLOSE(evalEnRidge, evalRidge, eps);
  BOOST_CHECK_CLOSE(evalEnLasso, evalLasso, eps);
  BOOST_CHECK_CLOSE(evalEnZero, evalZero, eps);
  BOOST_CHECK_CLOSE(evalEnMixed, evalMixed, eps);

  for (size_t i = 0; i < testVec.getSize(); ++i) {
    BOOST_CHECK_CLOSE(proxEnRidge[i], proxRidge[0], eps);
    BOOST_CHECK_CLOSE(proxEnLasso[i], proxLasso[0], eps);
    BOOST_CHECK_CLOSE(proxEnZero[i], proxZero[0], eps);
    BOOST_CHECK_CLOSE(proxEnMixed[i], proxMixed[0], eps);
  }
}

BOOST_AUTO_TEST_CASE(testLassoFunction) {
  auto fun = LassoFunction(lambda);
  const double evalIs = fun.eval(testVec);
  const auto proxIs = fun.prox(testVec, stepsize);

  const double evalShould = 10;
  const auto proxShould = DataVector(5, -1.0);

  BOOST_CHECK_CLOSE(evalIs, evalShould, eps);
  for (size_t i = 0; i < testVec.getSize(); ++i) {
    BOOST_CHECK_CLOSE(proxIs[i], proxShould[i], eps);
  }
}

BOOST_AUTO_TEST_CASE(testRidgeFunction) {
  auto fun = RidgeFunction(lambda);
  const double evalIs = fun.eval(testVec);
  const auto proxIs = fun.prox(testVec, stepsize);

  const double evalShould = 20;
  const auto proxShould = DataVector(5, -2.0 / 3.0);

  BOOST_CHECK_CLOSE(evalIs, evalShould, eps);
  for (size_t i = 0; i < testVec.getSize(); ++i) {
    BOOST_CHECK_CLOSE(proxIs[i], proxShould[i], eps);
  }
}

BOOST_AUTO_TEST_SUITE_END()
