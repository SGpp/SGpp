// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE SGppQuadratureModule
#include <boost/test/unit_test.hpp>

// fix for clang (from https://stackoverflow.com/a/33755176)
#ifdef __clang__
#include <string>

namespace boost {
namespace unit_test {
namespace ut_detail {

std::string normalize_test_case_name(const_string name) {
  return ((name[0] == '&') ? std::string(name.begin() + 1, name.size() - 1)
                           : std::string(name.begin(), name.size()));
}

}  // namespace ut_detail
}  // namespace unit_test
}  // namespace boost
#endif

#include <sgpp_base.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>

#include <sgpp_quadrature.hpp>
#include <sgpp/quadrature/QuadratureOpFactory.hpp>
#include <sgpp/globaldef.hpp>

#include <vector>

using sgpp::base::DataVector;
using sgpp::base::Grid;
using sgpp::quadrature::HaltonSampleGenerator;
using sgpp::quadrature::LatinHypercubeSampleGenerator;
using sgpp::quadrature::NaiveSampleGenerator;
using sgpp::quadrature::SampleGenerator;
using sgpp::quadrature::StratifiedSampleGenerator;

double f(DataVector x) {
  double res = 1.0f;

  for (size_t i = 0; i < x.getSize(); i++) {
    res *= 4 * (1 - x[i]) * x[i];
  }

  return res;
}

void testSampler(SampleGenerator& sampler, size_t dim, size_t numSamples, double analyticResult,
                 double tol) {
  DataVector sample = DataVector(dim);
  double sum = 0.0f;

  for (size_t i = 0; i < numSamples; i++) {
    sampler.getSample(sample);
    sum += f(sample);
  }

  sum /= static_cast<double>(numSamples);
  BOOST_CHECK_CLOSE(sum, analyticResult, tol * 1e2);
}

BOOST_AUTO_TEST_CASE(testSamplers) {
  size_t dim = 2;
  size_t numSamples = 100000;
  double analyticResult = std::pow(2. / 3., dim);
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
                                       sgpp::quadrature::SamplerTypes samplerType, size_t dim,
                                       size_t numSamples, std::vector<size_t>& blockSize,
                                       double analyticResult, double tol, uint64_t seed = 1234567) {
  std::unique_ptr<sgpp::quadrature::OperationQuadratureMCAdvanced> opQuad(
      sgpp::op_factory::createOperationQuadratureMCAdvanced(grid, numSamples, seed));

  switch (samplerType) {
    case sgpp::quadrature::SamplerTypes::Naive:
      opQuad->useNaiveMonteCarlo();
      break;

    case sgpp::quadrature::SamplerTypes::Stratified:
      opQuad->useStratifiedMonteCarlo(blockSize);
      break;

    case sgpp::quadrature::SamplerTypes::LatinHypercube:
      opQuad->useLatinHypercubeMonteCarlo();
      break;

    case sgpp::quadrature::SamplerTypes::Halton:
      opQuad->useQuasiMonteCarloWithHaltonSequences();
      break;

    default:
      std::cout << "test_quadrature::testOperationQuadratureMCAdvanced : sampler type not available"
                << std::endl;
  }

  double resMC = opQuad->doQuadrature(alpha);
  BOOST_CHECK_CLOSE(resMC, analyticResult, tol * 1e2);
}

BOOST_AUTO_TEST_CASE(testOperationMCAdvanced) {
  size_t dim = 2;
  size_t numSamples = 100000;
  double analyticResult = std::pow(2. / 3., dim);
  std::uint64_t seed = 1234567;

  // interpolate f on a sparse grid
  std::unique_ptr<sgpp::base::Grid> grid(sgpp::base::Grid::createPolyGrid(dim, 2));
  grid->getGenerator().regular(1);

  DataVector alpha(1);
  alpha[0] = 1.0f;

  std::unique_ptr<sgpp::base::OperationQuadrature> opQuad(
      sgpp::op_factory::createOperationQuadrature(*grid));
  double analyticIntegral = opQuad->doQuadrature(alpha);

  BOOST_CHECK_CLOSE(analyticIntegral, analyticResult, 1e-12);

  std::vector<size_t> blockSize(dim);

  for (size_t i = 0; i < dim; i++) {
    blockSize[i] = 10;
  }

  testOperationQuadratureMCAdvanced(*grid, alpha, sgpp::quadrature::SamplerTypes::Naive, dim,
                                    numSamples, blockSize, analyticResult, 5e-2, seed);
  testOperationQuadratureMCAdvanced(*grid, alpha, sgpp::quadrature::SamplerTypes::Stratified, dim,
                                    numSamples, blockSize, analyticResult, 1e-3, seed);
  testOperationQuadratureMCAdvanced(*grid, alpha, sgpp::quadrature::SamplerTypes::LatinHypercube,
                                    dim, numSamples, blockSize, analyticResult, 1e-3, seed);
  testOperationQuadratureMCAdvanced(*grid, alpha, sgpp::quadrature::SamplerTypes::Halton, dim,
                                    numSamples, blockSize, analyticResult, 1e-3, seed);
}
