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

#include <sgpp/combigrid/definitions.hpp>
#include <sgpp/globaldef.hpp>

#include <vector>

namespace ishigami_params {
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
}  // namespace ishigami_params

double ishigami_function(sgpp::base::DataVector const &v) {
  // transform [0, 1] -> [-pi, pi]
  sgpp::base::DataVector x(v);
  x.mult(2 * M_PI);
  x.sub(ishigami_params::pi);

  // evaluate the Ishigami function
  return std::sin(x[0]) + ishigami_params::a * std::sin(x[1]) * std::sin(x[1]) +
         ishigami_params::b * x[2] * x[2] * x[2] * x[2] * std::sin(x[0]);
}

#ifdef USE_DAKOTA
BOOST_AUTO_TEST_SUITE(testPolynomialChaosExpansion)

void testPCEMomentsAndSobolIndices(
    std::shared_ptr<sgpp::combigrid::CombigridOperation> op,
    std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D> functionBasis) {
  // compute variance of the estimator
  sgpp::combigrid::PolynomialChaosExpansion pce(op, functionBasis);

  // check the moments
  BOOST_CHECK_SMALL(std::abs(ishigami_params::mean - pce.mean()), ishigami_params::tolerance);
  BOOST_CHECK_SMALL(std::abs(ishigami_params::variance - pce.variance()),
                    ishigami_params::tolerance);

  // check the sobol indices
  sgpp::base::DataVector sobolIndices;
  pce.getComponentSobolIndices(sobolIndices);
  for (size_t i = 0; i < sobolIndices.size(); i++) {
    BOOST_CHECK_SMALL(std::abs(ishigami_params::sobolIndices[i] - sobolIndices[i]),
                      ishigami_params::tolerance);
  }

  // check the total sobol indices
  sgpp::base::DataVector totalSobolIndices;
  pce.getTotalSobolIndices(totalSobolIndices);
  for (size_t i = 0; i < totalSobolIndices.size(); i++) {
    BOOST_CHECK_SMALL(std::abs(ishigami_params::totalSobolIndices[i] - totalSobolIndices[i]),
                      ishigami_params::tolerance);
  }
}

BOOST_AUTO_TEST_CASE(testPCE) {
  sgpp::combigrid::OrthogonalPolynomialBasis1DConfiguration config;
  config.polyParameters.type_ = sgpp::combigrid::OrthogonalPolynomialBasisType::LEGENDRE;
  auto functionBasis = std::make_shared<sgpp::combigrid::OrthogonalPolynomialBasis1D>(config);

  sgpp::combigrid::MultiFunction func(ishigami_function);
  size_t level = 5;
  auto op = sgpp::combigrid::CombigridOperation::createExpL2LejaPolynomialInterpolation(
      ishigami_params::numDims, func);
  op->getLevelManager()->addRegularLevels(level);
  testPCEMomentsAndSobolIndices(op, functionBasis);
}

#endif

BOOST_AUTO_TEST_SUITE_END()
