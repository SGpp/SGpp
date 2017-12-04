// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp/combigrid/operation/onedim/PolynomialQuadratureEvaluator.hpp>
#include <sgpp/combigrid/operation/CombigridTensorOperation.hpp>
#include <sgpp/combigrid/functions/OrthogonalPolynomialBasis1D.hpp>
#include <sgpp/quadrature/sampling/LatinHypercubeSampleGenerator.hpp>
#include <sgpp/combigrid/pce/PolynomialChaosExpansion.hpp>
#include <sgpp/combigrid/pce/PolynomialStochasticCollocation.hpp>
#include <sgpp/combigrid/utils/AnalyticModels.hpp>
#include <sgpp/combigrid/definitions.hpp>
#include <sgpp/globaldef.hpp>

#include <vector>

#ifdef USE_DAKOTA

BOOST_AUTO_TEST_SUITE(testPolynomialChaosExpansion)

void testPCEIshigami(std::shared_ptr<sgpp::combigrid::CombigridOperation> op,
                     std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D> functionBasis) {
  // compute variance of the estimator
  sgpp::combigrid::PolynomialChaosExpansion pce(op, functionBasis);

  // ishigami model
  sgpp::combigrid::Ishigami ishigamiModel;

  // check the moments
  BOOST_CHECK_SMALL(std::abs(ishigamiModel.mean - pce.mean()), ishigamiModel.tolerance);
  BOOST_CHECK_SMALL(std::abs(ishigamiModel.variance - pce.variance()), ishigamiModel.tolerance);

  // check the sobol indices
  sgpp::base::DataVector sobolIndices;
  pce.getComponentSobolIndices(sobolIndices);
  for (size_t i = 0; i < sobolIndices.getSize(); i++) {
    BOOST_CHECK_SMALL(std::abs(ishigamiModel.sobolIndices[i] - sobolIndices[i]),
                      ishigamiModel.tolerance);
  }

  // check the total sobol indices
  sgpp::base::DataVector totalSobolIndices;
  pce.getTotalSobolIndices(totalSobolIndices);
  for (size_t i = 0; i < totalSobolIndices.getSize(); i++) {
    BOOST_CHECK_SMALL(std::abs(ishigamiModel.totalSobolIndices[i] - totalSobolIndices[i]),
                      ishigamiModel.tolerance);
  }
}

BOOST_AUTO_TEST_CASE(testPCE) {
  sgpp::combigrid::OrthogonalPolynomialBasis1DConfiguration config;
  config.polyParameters.type_ = sgpp::combigrid::OrthogonalPolynomialBasisType::LEGENDRE;
  auto functionBasis = std::make_shared<sgpp::combigrid::OrthogonalPolynomialBasis1D>(config);

  sgpp::combigrid::Ishigami ishigamiModel;
  sgpp::combigrid::MultiFunction func(ishigamiModel.eval);
  size_t level = 5;
  auto op = sgpp::combigrid::CombigridOperation::createExpL2LejaPolynomialInterpolation(
      ishigamiModel.numDims, func);
  op->getLevelManager()->addRegularLevels(level);
  testPCEIshigami(op, functionBasis);
}

void testPCEParbola(
    std::shared_ptr<sgpp::combigrid::CombigridOperation> op,
    std::vector<std::shared_ptr<sgpp::combigrid::AbstractInfiniteFunctionBasis1D>>& functionBases) {
  // compute variance of the estimator
  sgpp::combigrid::PolynomialChaosExpansion pce(op, functionBases);

  // check the moments
  sgpp::combigrid::Parabola parabolaModel;
  BOOST_CHECK_SMALL(std::abs(parabolaModel.mean - pce.mean()), parabolaModel.tolerance);
  BOOST_CHECK_SMALL(std::abs(parabolaModel.variance - pce.variance()), parabolaModel.tolerance);
}

BOOST_AUTO_TEST_CASE(testPCE_parabola) {
  sgpp::combigrid::Parabola parabolaModel;

  sgpp::combigrid::OrthogonalPolynomialBasis1DConfiguration config;
  config.polyParameters.type_ = sgpp::combigrid::OrthogonalPolynomialBasisType::JACOBI;
  config.polyParameters.lowerBound_ = parabolaModel.bounds[0];
  config.polyParameters.upperBound_ = parabolaModel.bounds[1];
  config.polyParameters.alpha_ = parabolaModel.alpha1;
  config.polyParameters.beta_ = parabolaModel.beta1;
  auto functionBasis1 = std::make_shared<sgpp::combigrid::OrthogonalPolynomialBasis1D>(config);

  config.polyParameters.type_ = sgpp::combigrid::OrthogonalPolynomialBasisType::JACOBI;
  config.polyParameters.lowerBound_ = parabolaModel.bounds[0];
  config.polyParameters.upperBound_ = parabolaModel.bounds[1];
  config.polyParameters.alpha_ = parabolaModel.alpha2;
  config.polyParameters.beta_ = parabolaModel.beta2;
  auto functionBasis2 = std::make_shared<sgpp::combigrid::OrthogonalPolynomialBasis1D>(config);

  std::vector<std::shared_ptr<sgpp::combigrid::AbstractInfiniteFunctionBasis1D>> functionBases{
      functionBasis1, functionBasis2};

  sgpp::combigrid::MultiFunction func(parabolaModel.eval);
  size_t level = 2;
  auto op = sgpp::combigrid::CombigridOperation::createExpL2LejaPolynomialInterpolation(
      parabolaModel.numDims, func);
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
  sgpp::combigrid::Ishigami ishigamiModel;
  sgpp::combigrid::OrthogonalPolynomialBasis1DConfiguration config;
  config.polyParameters.type_ = sgpp::combigrid::OrthogonalPolynomialBasisType::LEGENDRE;
  config.polyParameters.lowerBound_ = ishigamiModel.bounds[0];
  config.polyParameters.upperBound_ = ishigamiModel.bounds[1];
  auto functionBasis = std::make_shared<sgpp::combigrid::OrthogonalPolynomialBasis1D>(config);

  sgpp::combigrid::MultiFunction func(ishigamiModel.eval);
  size_t level = 5;
  auto op = sgpp::combigrid::CombigridOperation::createExpL2LejaPolynomialInterpolation(
      ishigamiModel.numDims, func);
  op->getLevelManager()->addRegularLevels(level);
  testStochasticCollocationMoments_id_marginals(op, functionBasis, ishigamiModel.mean,
                                                ishigamiModel.variance, ishigamiModel.tolerance);
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
  sgpp::combigrid::Parabola parabolaModel;
  sgpp::combigrid::OrthogonalPolynomialBasis1DConfiguration config;
  config.polyParameters.type_ = sgpp::combigrid::OrthogonalPolynomialBasisType::JACOBI;
  config.polyParameters.lowerBound_ = parabolaModel.bounds[0];
  config.polyParameters.upperBound_ = parabolaModel.bounds[1];
  config.polyParameters.alpha_ = parabolaModel.alpha1;
  config.polyParameters.beta_ = parabolaModel.beta1;
  auto functionBasis1 = std::make_shared<sgpp::combigrid::OrthogonalPolynomialBasis1D>(config);

  config.polyParameters.type_ = sgpp::combigrid::OrthogonalPolynomialBasisType::JACOBI;
  config.polyParameters.lowerBound_ = parabolaModel.bounds[0];
  config.polyParameters.upperBound_ = parabolaModel.bounds[1];
  config.polyParameters.alpha_ = parabolaModel.alpha2;
  config.polyParameters.beta_ = parabolaModel.beta2;
  auto functionBasis2 = std::make_shared<sgpp::combigrid::OrthogonalPolynomialBasis1D>(config);

  std::vector<std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D>> functionBases{
      functionBasis1, functionBasis2};

  sgpp::combigrid::MultiFunction func(parabolaModel.eval);
  size_t level = 2;
  auto op = sgpp::combigrid::CombigridOperation::createExpL2LejaPolynomialInterpolation(
      parabolaModel.numDims, func);
  op->getLevelManager()->addRegularLevels(level);
  testStochasticCollocationMoments_various_marginals(
      op, functionBases, parabolaModel.mean, parabolaModel.variance, parabolaModel.tolerance);
}

BOOST_AUTO_TEST_CASE(testStochasticCollocation_co2_lognormal) {
  sgpp::combigrid::CO2 co2Model;
  sgpp::combigrid::OrthogonalPolynomialBasis1DConfiguration config;
  config.polyParameters.type_ = sgpp::combigrid::OrthogonalPolynomialBasisType::BOUNDED_LOGNORMAL;
  config.polyParameters.logmean_ = co2Model.logmean;
  config.polyParameters.stddev_ = co2Model.stddev;
  config.polyParameters.lowerBound_ = co2Model.bounds[0];
  config.polyParameters.upperBound_ = co2Model.bounds[1];
  auto functionBasis = std::make_shared<sgpp::combigrid::OrthogonalPolynomialBasis1D>(config);
  sgpp::combigrid::MultiFunction func(co2Model.eval);

  auto op = sgpp::combigrid::CombigridOperation::createExpClenshawCurtisPolynomialInterpolation(
      co2Model.numDims, func);
  auto op_levelManager = op->getLevelManager();

  // compute variance of the estimator
  sgpp::combigrid::PolynomialStochasticCollocation sc(op, functionBasis);
  auto tensor_levelManager = sc.getCombigridTensorOperation()->getLevelManager();

  op_levelManager->addRegularLevels(2);
  tensor_levelManager->addLevelsFromStructure(op_levelManager->getLevelStructure());

  // check the moments
  BOOST_CHECK_SMALL(std::abs(co2Model.mean - sc.mean()), 1e-1);
  BOOST_CHECK_SMALL(std::abs(co2Model.variance - sc.variance()), 1e-1);

  auto op2 = sgpp::combigrid::CombigridOperation::createExpL2LejaPolynomialInterpolation(
      co2Model.numDims, func);
  auto op2_levelManager = op2->getLevelManager();

  // compute variance of the estimator
  sc.updateOperation(op2);
  tensor_levelManager = sc.getCombigridTensorOperation()->getLevelManager();

  op2_levelManager->addRegularLevels(7);
  tensor_levelManager->addLevelsFromStructure(op2_levelManager->getLevelStructure());

  // check the moments
  BOOST_CHECK_SMALL(std::abs(co2Model.mean - sc.mean()), co2Model.tolerance);
  BOOST_CHECK_SMALL(std::abs(co2Model.variance - sc.variance()), co2Model.tolerance);
}

BOOST_AUTO_TEST_SUITE_END()

#endif
