// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE SGppQuadratureModule
#include <boost/test/unit_test.hpp>

#include <sgpp_base.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>

#include <sgpp_quadrature.hpp>
#include <sgpp/quadrature/QuadratureOpFactory.hpp>
#include <sgpp/globaldef.hpp>

#include <vector>

using SGPP::base::DataVector;
using SGPP::base::Grid;
using SGPP::quadrature::HaltonSampleGenerator;
using SGPP::quadrature::LatinHypercubeSampleGenerator;
using SGPP::quadrature::NaiveSampleGenerator;
using SGPP::quadrature::SampleGenerator;
using SGPP::quadrature::StratifiedSampleGenerator;

SGPP::float_t f(DataVector x) {
  SGPP::float_t res = 1.0f;

  for (size_t i = 0; i < x.getSize(); i++) {
    res *= 4 * (1 - x[i]) * x[i];
  }

  return res;
}

void testSampler(SampleGenerator& sampler, size_t dim, size_t numSamples,
                 SGPP::float_t analyticResult, double tol) {
  DataVector sample = DataVector(dim);
  SGPP::float_t sum = 0.0f;

  for (size_t i = 0; i < numSamples; i++) {
    sampler.getSample(sample);
    sum += f(sample);
  }

  sum /= static_cast<SGPP::float_t>(numSamples);
  BOOST_CHECK_CLOSE(sum, analyticResult, tol * 1e2);
}

BOOST_AUTO_TEST_CASE(testSamplers) {
  size_t dim = 2;
  size_t numSamples = 100000;
  SGPP::float_t analyticResult = std::pow(2. / 3., dim);
  uint64_t seed = 1234567;

  NaiveSampleGenerator pNSampler(dim, seed);
  HaltonSampleGenerator pHSampler(dim);
  LatinHypercubeSampleGenerator pLHSampler(dim, numSamples, seed);
  std::vector<size_t> blockSize(dim, seed);

  for (size_t i = 0; i < dim; i++) {
    blockSize[i] = 10;
  }

  StratifiedSampleGenerator pSSampler(blockSize);

  testSampler(pNSampler, dim, numSamples, analyticResult, 5e-2);
  testSampler(pHSampler, dim, numSamples, analyticResult, 1e-3);
  testSampler(pLHSampler, dim, numSamples, analyticResult, 1e-3);
  testSampler(pSSampler, dim, numSamples, analyticResult, 1e-3);
}

void testOperationQuadratureMCAdvanced(Grid& grid, DataVector& alpha,
                                       SGPP::quadrature::SamplerTypes samplerType, size_t dim,
                                       size_t numSamples,
                                       std::vector<size_t>& blockSize, SGPP::float_t analyticResult,
                                       double tol, uint64_t seed = 1234567) {
  SGPP::quadrature::OperationQuadratureMCAdvanced* opQuad =
    SGPP::op_factory::createOperationQuadratureMCAdvanced(grid, numSamples, seed);

  switch (samplerType) {
    case SGPP::quadrature::SamplerTypes::Naive:
      opQuad->useNaiveMonteCarlo();
      break;

    case SGPP::quadrature::SamplerTypes::Stratified:
      opQuad->useStratifiedMonteCarlo(blockSize);
      break;

    case SGPP::quadrature::SamplerTypes::LatinHypercube:
      opQuad->useLatinHypercubeMonteCarlo();
      break;

    case SGPP::quadrature::SamplerTypes::Halton:
      opQuad->useQuasiMonteCarloWithHaltonSequences();
      break;

    default:
      std::cout
          << "test_quadrature::testOperationQuadratureMCAdvanced : sampler type not available"
          << std::endl;
  }

  SGPP::float_t resMC = opQuad->doQuadrature(alpha);
  BOOST_CHECK_CLOSE(resMC, analyticResult, tol * 1e2);
}

BOOST_AUTO_TEST_CASE(testOperationMCAdvanced) {
  size_t dim = 2;
  size_t numSamples = 100000;
  SGPP::float_t analyticResult = std::pow(2. / 3., dim);
  std::uint64_t seed = 1234567;

  // interpolate f on a sparse grid
  SGPP::base::Grid* grid = SGPP::base::Grid::createPolyGrid(dim, 2);
  SGPP::base::GridGenerator* gridGen = grid->createGridGenerator();
  gridGen->regular(1);

  DataVector alpha(1);
  alpha[0] = 1.0f;

  SGPP::base::OperationQuadrature* opQuad =
    SGPP::op_factory::createOperationQuadrature(*grid);
  SGPP::float_t analyticIntegral = opQuad->doQuadrature(alpha);

#if USE_DOUBLE_PRECISION == 1
  BOOST_CHECK_CLOSE(analyticIntegral, analyticResult, float_t(1e-12) );
#else
  BOOST_CHECK_CLOSE(analyticIntegral, analyticResult, float_t(1e-5) );
#endif

  std::vector<size_t> blockSize(dim);

  for (size_t i = 0; i < dim; i++) {
    blockSize[i] = 10;
  }

  testOperationQuadratureMCAdvanced(*grid, alpha,
                                    SGPP::quadrature::SamplerTypes::Naive, dim, numSamples,
                                    blockSize,
                                    analyticResult, 5e-2, seed);
  testOperationQuadratureMCAdvanced(*grid, alpha,
                                    SGPP::quadrature::SamplerTypes::Stratified, dim, numSamples,
                                    blockSize,
                                    analyticResult, 1e-3, seed);
  testOperationQuadratureMCAdvanced(*grid, alpha,
                                    SGPP::quadrature::SamplerTypes::LatinHypercube, dim, numSamples,
                                    blockSize,
                                    analyticResult, 1e-3, seed);
  testOperationQuadratureMCAdvanced(*grid, alpha,
                                    SGPP::quadrature::SamplerTypes::Halton, dim, numSamples,
                                    blockSize,
                                    analyticResult, 1e-3, seed);
}
