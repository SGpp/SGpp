// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp/combigrid/functions/OrthogonalPolynomialBasis1D.hpp>
#include <sgpp/combigrid/operation/CombigridTensorOperation.hpp>
#include <sgpp/combigrid/operation/onedim/PolynomialQuadratureEvaluator.hpp>
#include <sgpp/combigrid/pce/CombigridSurrogateModel.hpp>
#include <sgpp/combigrid/pce/CombigridSurrogateModelFactory.hpp>
#include <sgpp/combigrid/utils/AnalyticModels.hpp>
#include <sgpp/quadrature/sampling/LatinHypercubeSampleGenerator.hpp>

#include <vector>

const double tolerance = 1e-12;

BOOST_AUTO_TEST_SUITE(testIntegration)

BOOST_AUTO_TEST_CASE(testQuadrature) {
  sgpp::combigrid::PolynomialQuadratureEvaluator eval;

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
  sgpp::combigrid::PolynomialQuadratureEvaluator eval;

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
  eval.setBasisCoefficientsAtGridPoints(function_points);

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
  eval.setBasisCoefficientsAtGridPoints(function_points);

  BOOST_CHECK_CLOSE(eval.eval().getValue(), accurate_solution, 0.2);
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

double mean(size_t numDims, sgpp::combigrid::MultiFunction& func, size_t numPoints,
            std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D> functionBasis) {
  sgpp::combigrid::MultiFunction mean_func(
      [func, functionBasis](sgpp::base::DataVector const& param) {
        double value = func(param);
        double pdf_value = 1.0;
        for (size_t i = 0; i < param.getSize(); i++) {
          pdf_value *= functionBasis->pdf(param[i]);
        }
        return value * pdf_value;
      });
  return monte_carlo_quadrature(numDims, mean_func, numPoints);
}

double variance(size_t numDims, sgpp::combigrid::MultiFunction& func, size_t numPoints,
                double mean_ref,
                std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D> functionBasis) {
  sgpp::combigrid::MultiFunction var_func(
      [func, functionBasis, mean_ref](sgpp::base::DataVector const& param) {
        double value = std::pow(func(param) - mean_ref, 2);
        double pdf_value = 1.0;
        for (size_t i = 0; i < param.getSize(); i++) {
          pdf_value *= functionBasis->pdf(param[i]);
        }
        return value * pdf_value;
      });
  return monte_carlo_quadrature(numDims, var_func, numPoints);
}

#ifdef USE_DAKOTA
void testTensorOperation(
    std::shared_ptr<sgpp::combigrid::CombigridTensorOperation> tensor_op,
    std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D> basisFunction, size_t level,
    size_t numDims) {
  // compute variance of the estimator
  tensor_op->getLevelManager()->addRegularLevels(level);

  // initialize the surrogate model
  sgpp::combigrid::CombigridSurrogateModelConfiguration config;
  config.type = sgpp::combigrid::CombigridSurrogateModelsType::POLYNOMIAL_CHAOS_EXPANSION;
  config.loadFromCombigridOperation(tensor_op);
  config.basisFunction = basisFunction;
  auto pce = sgpp::combigrid::createCombigridSurrogateModel(config);

  // compute the reference solution
  double mean_ref = sgpp::combigrid::Parabola_uniform::mean(numDims);
  double var_ref = sgpp::combigrid::Parabola_uniform::variance(numDims);

  // check if they match
  BOOST_CHECK_CLOSE(pce->mean(), mean_ref, 1e-10);
  BOOST_CHECK_CLOSE(pce->variance(), var_ref, 1e-10);
}

BOOST_AUTO_TEST_CASE(testVarianceComputationWithPCETransformation_expLeja) {
  sgpp::combigrid::OrthogonalPolynomialBasis1DConfiguration config;
  config.polyParameters.type_ = sgpp::combigrid::OrthogonalPolynomialBasisType::LEGENDRE;
  auto functionBasis = std::make_shared<sgpp::combigrid::OrthogonalPolynomialBasis1D>(config);

  sgpp::combigrid::MultiFunction func(sgpp::combigrid::Parabola_uniform::eval);
  for (size_t numDims = 2; numDims <= 5; ++numDims) {
    auto tensor_op =
        sgpp::combigrid::CombigridTensorOperation::createExpLejaPolynomialInterpolation(
            functionBasis, numDims, func);
    testTensorOperation(tensor_op, functionBasis, numDims + 1, numDims);
  }
}

BOOST_AUTO_TEST_CASE(testVarianceComputationWithPCETransformation_expL2Leja) {
  sgpp::combigrid::OrthogonalPolynomialBasis1DConfiguration config;
  config.polyParameters.type_ = sgpp::combigrid::OrthogonalPolynomialBasisType::LEGENDRE;
  auto functionBasis = std::make_shared<sgpp::combigrid::OrthogonalPolynomialBasis1D>(config);

  sgpp::combigrid::MultiFunction func(sgpp::combigrid::Parabola_uniform::eval);
  for (size_t numDims = 2; numDims <= 5; ++numDims) {
    auto tensor_op =
        sgpp::combigrid::CombigridTensorOperation::createExpL2LejaPolynomialInterpolation(
            functionBasis, numDims, func);
    testTensorOperation(tensor_op, functionBasis, numDims + 1, numDims);
  }
}

BOOST_AUTO_TEST_CASE(testVarianceComputationWithPCETransformation_ClenshawCurtis) {
  sgpp::combigrid::OrthogonalPolynomialBasis1DConfiguration config;
  config.polyParameters.type_ = sgpp::combigrid::OrthogonalPolynomialBasisType::LEGENDRE;
  auto functionBasis = std::make_shared<sgpp::combigrid::OrthogonalPolynomialBasis1D>(config);

  sgpp::combigrid::MultiFunction func(sgpp::combigrid::Parabola_uniform::eval);
  for (size_t numDims = 2; numDims <= 5; ++numDims) {
    auto tensor_op =
        sgpp::combigrid::CombigridTensorOperation::createExpClenshawCurtisPolynomialInterpolation(
            functionBasis, numDims, func);
    testTensorOperation(tensor_op, functionBasis, numDims + 1, numDims);
  }
}

BOOST_AUTO_TEST_CASE(testVarianceComputationWithPCETransformation_Lognormal_ClenshawCurtis) {
  sgpp::combigrid::OrthogonalPolynomialBasis1DConfiguration config;
  config.polyParameters.type_ = sgpp::combigrid::OrthogonalPolynomialBasisType::BOUNDED_LOGNORMAL;
  config.polyParameters.logmean_ = 0.0;
  config.polyParameters.stddev_ = 1.0;
  config.polyParameters.lowerBound_ = 1e-2;
  config.polyParameters.upperBound_ = 1.0;
  auto functionBasis = std::make_shared<sgpp::combigrid::OrthogonalPolynomialBasis1D>(config);

  sgpp::combigrid::CombigridSurrogateModelConfiguration pce_config;
  pce_config.type = sgpp::combigrid::CombigridSurrogateModelsType::POLYNOMIAL_CHAOS_EXPANSION;
  pce_config.basisFunction = functionBasis;

  sgpp::combigrid::MultiFunction func(sgpp::combigrid::Parabola_uniform::eval);
  size_t level = 4;

  sgpp::combigrid::CombigridSurrogateModelConfiguration pce_update_config;
  for (size_t numDims = 1; numDims <= 5; ++numDims) {
    auto op = sgpp::combigrid::CombigridOperation::createExpClenshawCurtisPolynomialInterpolation(
        numDims, func);
    auto op_levelManager = op->getLevelManager();

    // initialize the surrogate model
    pce_config.loadFromCombigridOperation(op);
    auto pce = sgpp::combigrid::createCombigridSurrogateModel(pce_config);

    // compute the reference solution
    double mean_ref = mean(numDims, func, 1e4, functionBasis);
    double var_ref = variance(numDims, func, 1e4, mean_ref, functionBasis);
    op_levelManager->addRegularLevels(level);

    // update pce
    pce_update_config.levelStructure = op_levelManager->getLevelStructure();
    pce->updateConfig(pce_update_config);

    // check if they match
    //    std::cout << "  level = " << level << ": mean = " << mean_ref << " ~ " << pce.mean()
    //              << " (err = " << std::abs(pce.mean() - mean_ref) << "); var = " << var_ref << "
    //              ~ "
    //              << pce.variance() << " (err = " << std::abs(pce.variance() - var_ref) << ")"
    //              << std::endl;
    BOOST_CHECK_SMALL(std::abs(pce->mean() - mean_ref), 5e-3);
    BOOST_CHECK_SMALL(std::abs(pce->variance() - var_ref), 5e-3);
  }
}

#endif

BOOST_AUTO_TEST_SUITE_END()
