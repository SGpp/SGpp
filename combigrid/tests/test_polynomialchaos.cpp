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
#include <sgpp/combigrid/pce/PolynomialChaosExpansion.hpp>
#include <sgpp/combigrid/pce/PolynomialStochasticCollocation.hpp>

#include <sgpp/combigrid/definitions.hpp>
#include <sgpp/globaldef.hpp>

#include <vector>

// ----------------------------------------------------------------------------------
namespace ishigami {
double tolerance = 5e-4;
size_t numDims = 3;
std::vector<double> pi_vec{M_PI, M_PI, M_PI};
sgpp::base::DataVector pi(pi_vec);
double pi_4 = M_PI * M_PI * M_PI * M_PI;
double a = 7.0;
double b = 0.1;
double variance = a * a / 8. + b * pi_4 / 5 + b * b * pi_4 * pi_4 / 18. + 0.5;
double mean = 3.5;
std::vector<double> sobolIndices{0.3138, 0.4424, 0.0, 0.0, 0.2436, 0.0, 0.0};
std::vector<double> totalSobolIndices{0.5574, 0.4424, 0.2436};
std::vector<double> bounds{0, 1};

double eval(sgpp::base::DataVector const& v) {
  // transform [0, 1] -> [-pi, pi]
  sgpp::base::DataVector x(v);
  x.mult(2 * M_PI);
  x.sub(ishigami::pi);

  // evaluate the Ishigami function
  return std::sin(x[0]) + ishigami::a * std::sin(x[1]) * std::sin(x[1]) +
         ishigami::b * x[2] * x[2] * x[2] * x[2] * std::sin(x[0]);
}
}  // namespace ishigami
// ----------------------------------------------------------------------------------
namespace parabola {
double tolerance = 1e-13;
size_t numDims = 2;

double alpha1 = 5.0;
double beta1 = 4.0;

double alpha2 = 3.0;
double beta2 = 2.0;

double c1 = std::tgamma(alpha1 + beta1) / (std::tgamma(alpha1) * std::tgamma(beta1));
double c2 = std::tgamma(alpha2 + beta2) / (std::tgamma(alpha2) * std::tgamma(beta2));

double mean = c1 * c2 / 4725.0;
double variance = c1 * c1 * c1 * c2 * c2 * c2 / 75014100000.0 - 2 * c1 * c1 * c2 * c2 / 22325625.0 +
                  4.0 * c1 * c2 / 24255.0;
std::vector<double> bounds{0, 1};

double eval(sgpp::base::DataVector const& v) {
  double ans = 1.0;
  for (size_t idim = 0; idim < v.getSize(); idim++) {
    ans *= 4.0 * v[idim] * (1.0 - v[idim]);
  }
  return ans;
}
}  // namespace parabola
// ----------------------------------------------------------------------------------

#ifdef USE_DAKOTA

BOOST_AUTO_TEST_SUITE(testPolynomialChaosExpansion)

void testPCEIshigami(std::shared_ptr<sgpp::combigrid::CombigridOperation> op,
                     std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D> functionBasis) {
  // compute variance of the estimator
  sgpp::combigrid::PolynomialChaosExpansion pce(op, functionBasis);

  // check the moments
  BOOST_CHECK_SMALL(std::abs(ishigami::mean - pce.mean()), ishigami::tolerance);
  BOOST_CHECK_SMALL(std::abs(ishigami::variance - pce.variance()), ishigami::tolerance);

  // check the sobol indices
  sgpp::base::DataVector sobolIndices;
  pce.getComponentSobolIndices(sobolIndices);
  for (size_t i = 0; i < sobolIndices.size(); i++) {
    BOOST_CHECK_SMALL(std::abs(ishigami::sobolIndices[i] - sobolIndices[i]), ishigami::tolerance);
  }

  // check the total sobol indices
  sgpp::base::DataVector totalSobolIndices;
  pce.getTotalSobolIndices(totalSobolIndices);
  for (size_t i = 0; i < totalSobolIndices.size(); i++) {
    BOOST_CHECK_SMALL(std::abs(ishigami::totalSobolIndices[i] - totalSobolIndices[i]),
                      ishigami::tolerance);
  }
}

BOOST_AUTO_TEST_CASE(testPCE) {
  sgpp::combigrid::OrthogonalPolynomialBasis1DConfiguration config;
  config.polyParameters.type_ = sgpp::combigrid::OrthogonalPolynomialBasisType::LEGENDRE;
  auto functionBasis = std::make_shared<sgpp::combigrid::OrthogonalPolynomialBasis1D>(config);

  sgpp::combigrid::MultiFunction func(ishigami::eval);
  size_t level = 5;
  auto op = sgpp::combigrid::CombigridOperation::createExpL2LejaPolynomialInterpolation(
      ishigami::numDims, func);
  op->getLevelManager()->addRegularLevels(level);
  testPCEIshigami(op, functionBasis);
}

void testPCEParbola(
    std::shared_ptr<sgpp::combigrid::CombigridOperation> op,
    std::vector<std::shared_ptr<sgpp::combigrid::AbstractInfiniteFunctionBasis1D>>& functionBases) {
  // compute variance of the estimator
  sgpp::combigrid::PolynomialChaosExpansion pce(op, functionBases);

  // check the moments
  BOOST_CHECK_SMALL(std::abs(parabola::mean - pce.mean()), parabola::tolerance);
  BOOST_CHECK_SMALL(std::abs(parabola::variance - pce.variance()), parabola::tolerance);
}

BOOST_AUTO_TEST_CASE(testPCE_parabola) {
  sgpp::combigrid::OrthogonalPolynomialBasis1DConfiguration config;
  config.polyParameters.type_ = sgpp::combigrid::OrthogonalPolynomialBasisType::JACOBI;
  config.polyParameters.lowerBound_ = parabola::bounds[0];
  config.polyParameters.upperBound_ = parabola::bounds[1];
  config.polyParameters.alpha_ = parabola::alpha1;
  config.polyParameters.beta_ = parabola::beta1;
  auto functionBasis1 = std::make_shared<sgpp::combigrid::OrthogonalPolynomialBasis1D>(config);

  config.polyParameters.type_ = sgpp::combigrid::OrthogonalPolynomialBasisType::JACOBI;
  config.polyParameters.lowerBound_ = parabola::bounds[0];
  config.polyParameters.upperBound_ = parabola::bounds[1];
  config.polyParameters.alpha_ = parabola::alpha2;
  config.polyParameters.beta_ = parabola::beta2;
  auto functionBasis2 = std::make_shared<sgpp::combigrid::OrthogonalPolynomialBasis1D>(config);

  std::vector<std::shared_ptr<sgpp::combigrid::AbstractInfiniteFunctionBasis1D>> functionBases{
      functionBasis1, functionBasis2};

  sgpp::combigrid::MultiFunction func(parabola::eval);
  size_t level = 2;
  auto op = sgpp::combigrid::CombigridOperation::createExpL2LejaPolynomialInterpolation(
      parabola::numDims, func);
  op->getLevelManager()->addRegularLevels(level);
  testPCEParbola(op, functionBases);
}

BOOST_AUTO_TEST_SUITE_END()

// ----------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE(testPolynomialStochasticCollocation)

void testStochasticCollocationMoments_id_marginals(
    std::shared_ptr<sgpp::combigrid::CombigridOperation> op,
    std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D> functionBasis, double mean,
    double variance, double tol) {
  // compute variance of the estimator
  sgpp::combigrid::PolynomialStochasticCollocation sc(op, functionBasis);

  // check the moments
  BOOST_CHECK_SMALL(std::abs(mean - sc.mean()), tol);
  BOOST_CHECK_SMALL(std::abs(variance - sc.variance()), tol);
}

BOOST_AUTO_TEST_CASE(testStochasticCollocation_ishigami) {
  sgpp::combigrid::OrthogonalPolynomialBasis1DConfiguration config;
  config.polyParameters.type_ = sgpp::combigrid::OrthogonalPolynomialBasisType::LEGENDRE;
  config.polyParameters.lowerBound_ = ishigami::bounds[0];
  config.polyParameters.upperBound_ = ishigami::bounds[1];
  auto functionBasis = std::make_shared<sgpp::combigrid::OrthogonalPolynomialBasis1D>(config);

  sgpp::combigrid::MultiFunction func(ishigami::eval);
  size_t level = 5;
  auto op = sgpp::combigrid::CombigridOperation::createExpL2LejaPolynomialInterpolation(
      ishigami::numDims, func);
  op->getLevelManager()->addRegularLevels(level);
  testStochasticCollocationMoments_id_marginals(op, functionBasis, ishigami::mean,
                                                ishigami::variance, ishigami::tolerance);
}

void testStochasticCollocationMoments_various_marginals(
    std::shared_ptr<sgpp::combigrid::CombigridOperation> op,
    std::vector<std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D>>& functionBases,
    double mean, double variance, double tol) {
  // compute variance of the estimator
  sgpp::combigrid::PolynomialStochasticCollocation sc(op, functionBases);

  // check the moments
  BOOST_CHECK_SMALL(std::abs(mean - sc.mean()), tol);
  BOOST_CHECK_SMALL(std::abs(variance - sc.variance()), tol);
}

BOOST_AUTO_TEST_CASE(testStochasticCollocation_parabola) {
  sgpp::combigrid::OrthogonalPolynomialBasis1DConfiguration config;
  config.polyParameters.type_ = sgpp::combigrid::OrthogonalPolynomialBasisType::JACOBI;
  config.polyParameters.lowerBound_ = parabola::bounds[0];
  config.polyParameters.upperBound_ = parabola::bounds[1];
  config.polyParameters.alpha_ = parabola::alpha1;
  config.polyParameters.beta_ = parabola::beta1;
  auto functionBasis1 = std::make_shared<sgpp::combigrid::OrthogonalPolynomialBasis1D>(config);

  config.polyParameters.type_ = sgpp::combigrid::OrthogonalPolynomialBasisType::JACOBI;
  config.polyParameters.lowerBound_ = parabola::bounds[0];
  config.polyParameters.upperBound_ = parabola::bounds[1];
  config.polyParameters.alpha_ = parabola::alpha2;
  config.polyParameters.beta_ = parabola::beta2;
  auto functionBasis2 = std::make_shared<sgpp::combigrid::OrthogonalPolynomialBasis1D>(config);

  std::vector<std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D>> functionBases{
      functionBasis1, functionBasis2};

  sgpp::combigrid::MultiFunction func(parabola::eval);
  size_t level = 2;
  auto op = sgpp::combigrid::CombigridOperation::createExpL2LejaPolynomialInterpolation(
      parabola::numDims, func);
  op->getLevelManager()->addRegularLevels(level);
  testStochasticCollocationMoments_various_marginals(op, functionBases, parabola::mean,
                                                     parabola::variance, parabola::tolerance);
}

BOOST_AUTO_TEST_SUITE_END()

#endif
