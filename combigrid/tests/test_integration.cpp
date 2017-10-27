// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp/combigrid/operation/onedim/QuadratureEvaluator.hpp>
#include <sgpp/combigrid/operation/CombigridTensorOperation.hpp>
#include <sgpp/combigrid/functions/OrthogonalPolynomialBasis1D.hpp>
#include <sgpp/quadrature/sampling/LatinHypercubeSampleGenerator.hpp>

#include <sgpp/combigrid/definitions.hpp>
#include <sgpp/globaldef.hpp>

#include <vector>

const double tolerance = 1e-12;

BOOST_AUTO_TEST_SUITE(testIntegration)

BOOST_AUTO_TEST_CASE(testQuadrature) {
  sgpp::combigrid::QuadratureEvaluator eval;

  std::vector<double> points;
  std::vector<double> weights;

  // 1 Point

  points.push_back(0.0);

  eval.setGridPoints(points);

  BOOST_CHECK_EQUAL(eval.getBasisValues().size(), 1);
  BOOST_CHECK_CLOSE(eval.getBasisValues()[0].getValue(), 1.0, tolerance);

  // 2 Points

  points.push_back(1.0);

  eval.setGridPoints(points);

  BOOST_CHECK_EQUAL(eval.getBasisValues().size(), 2);

  BOOST_CHECK_CLOSE(eval.getBasisValues()[0].getValue(), 0.5, tolerance);
  BOOST_CHECK_CLOSE(eval.getBasisValues()[1].getValue(), 0.5, tolerance);
}

BOOST_AUTO_TEST_CASE(testQuadratureWithPolynomial) {
  sgpp::combigrid::QuadratureEvaluator eval;

  std::vector<double> points;
  std::vector<double> weights;

  points.push_back(0.0);
  points.push_back(0.5);
  points.push_back(1.0);

  eval.setGridPoints(points);

  // p(x) = x^2

  std::vector<double> function_points;
  function_points.push_back(0.0);
  function_points.push_back(0.25);
  function_points.push_back(1.0);
  eval.setFunctionValuesAtGridPoints(function_points);

  // check coefficients
  BOOST_CHECK_EQUAL(eval.getBasisValues().size(), 3);
  BOOST_CHECK_CLOSE(eval.getBasisValues()[0].getValue(), 1.0 / 6, tolerance);
  BOOST_CHECK_CLOSE(eval.getBasisValues()[1].getValue(), 2.0 / 3, tolerance);
  BOOST_CHECK_CLOSE(eval.getBasisValues()[2].getValue(), 1.0 / 6, tolerance);

  double accurate_solution = 1.0 / 3;

  BOOST_CHECK_CLOSE(eval.eval().getValue(), accurate_solution, tolerance);

  // try a more complex polynomial
  // p(x) = 6 * x^5 - 5 * x^4 + 4 * x^3 - 3 * x^2 + 2 * x - 42
  // the accurate integral is x * (x^5 - x^4 + x^3 - x^2 + x - 42)
  // so in [0,1], it is -41
  accurate_solution = -41.0;

  // we need to give d + 1 points, here d is 5
  points.clear();
  function_points.clear();
  points.push_back(0.0);
  function_points.push_back(-42.0);
  points.push_back(0.2);
  function_points.push_back(-41.69408);
  points.push_back(0.4);
  function_points.push_back(-41.49056);
  points.push_back(0.6);
  function_points.push_back(-41.19744);
  points.push_back(0.8);
  function_points.push_back(-40.35392);
  points.push_back(1.0);
  function_points.push_back(-38.0);

  eval.setGridPoints(points);
  eval.setFunctionValuesAtGridPoints(function_points);

  BOOST_CHECK_CLOSE(eval.eval().getValue(), accurate_solution, 0.2);
}

double normal_parabola(sgpp::base::DataVector const& x) {
  double ans = 1.0;
  for (size_t idim = 0; idim < x.getSize(); idim++) {
    ans *= 4 * x[idim] * (1.0 - x[idim]);
  }
  return ans;
}

double normal_parabola_mean_uniform(size_t numDims) { return std::pow(2. / 3., numDims); }

double normal_parabola_variance_uniform(size_t numDims) {
  double mean = normal_parabola_mean_uniform(numDims);
  return std::pow(16. / 30., numDims) - std::pow(mean, 2);
}

double monte_carlo_quadrature(size_t numDims, sgpp::combigrid::MultiFunction& func,
                              size_t numPoints = 10000) {
  sgpp::quadrature::LatinHypercubeSampleGenerator generator(numDims, numPoints);
  sgpp::base::DataVector p(numDims, 0);

  double mysum = 0.0;
  for (size_t i = 0; i < numPoints; i++) {
    generator.getSample(p);
    mysum += func(p);
  }

  return mysum / static_cast<double>(numPoints);
}

void testTensorOperation(std::shared_ptr<sgpp::combigrid::CombigridTensorOperation> tensor_op,
                         size_t level, size_t numDims) {
  // compute variance of the estimator
  auto tensor_result = tensor_op->evaluate(level);
  double mean = tensor_result.get(sgpp::combigrid::MultiIndex(numDims, 0)).getValue();
  double var = std::pow(tensor_result.norm(), 2);

  // compute the reference solution
  double mean_ref = normal_parabola_mean_uniform(numDims);
  double var_ref = normal_parabola_variance_uniform(numDims);

  // check if they match
  BOOST_CHECK_CLOSE(mean, mean_ref, 1e-10);
  BOOST_CHECK_CLOSE(var, var_ref, 1e-10);
}

BOOST_AUTO_TEST_CASE(testVarianceComputationWithPCETransformation_expLeja) {
  sgpp::combigrid::OrthogonalPolynomialBasis1DConfiguration config;
  config.polyParameters.type_ = sgpp::combigrid::OrthogonalPolynomialBasisType::LEGENDRE;
  auto functionBasis = std::make_shared<sgpp::combigrid::OrthogonalPolynomialBasis1D>(config);

  sgpp::combigrid::MultiFunction func(normal_parabola);
  for (size_t numDims = 2; numDims <= 5; ++numDims) {
    auto tensor_op =
        sgpp::combigrid::CombigridTensorOperation::createExpLejaPolynomialInterpolation(
            functionBasis, numDims, func);
    testTensorOperation(tensor_op, numDims + 1, numDims);
  }
}

BOOST_AUTO_TEST_CASE(testVarianceComputationWithPCETransformation_expL2Leja) {
  sgpp::combigrid::OrthogonalPolynomialBasis1DConfiguration config;
  config.polyParameters.type_ = sgpp::combigrid::OrthogonalPolynomialBasisType::LEGENDRE;
  auto functionBasis = std::make_shared<sgpp::combigrid::OrthogonalPolynomialBasis1D>(config);

  sgpp::combigrid::MultiFunction func(normal_parabola);
  for (size_t numDims = 2; numDims <= 5; ++numDims) {
    auto tensor_op =
        sgpp::combigrid::CombigridTensorOperation::createExpL2LejaPolynomialInterpolation(
            functionBasis, numDims, func);
    testTensorOperation(tensor_op, numDims + 1, numDims);
  }
}

BOOST_AUTO_TEST_CASE(testVarianceComputationWithPCETransformation_ClenshawCurtis) {
  sgpp::combigrid::OrthogonalPolynomialBasis1DConfiguration config;
  config.polyParameters.type_ = sgpp::combigrid::OrthogonalPolynomialBasisType::LEGENDRE;
  auto functionBasis = std::make_shared<sgpp::combigrid::OrthogonalPolynomialBasis1D>(config);

  sgpp::combigrid::MultiFunction func(normal_parabola);
  for (size_t numDims = 2; numDims <= 5; ++numDims) {
    auto tensor_op =
        sgpp::combigrid::CombigridTensorOperation::createExpClenshawCurtisPolynomialInterpolation(
            functionBasis, numDims, func);
    testTensorOperation(tensor_op, numDims + 1, numDims);
  }
}

BOOST_AUTO_TEST_SUITE_END()
