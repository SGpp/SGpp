// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp/combigrid/operation/onedim/PolynomialQuadratureEvaluator.hpp>
#include <sgpp/combigrid/operation/CombigridTensorOperation.hpp>
#include <sgpp/combigrid/functions/OrthogonalPolynomialBasis1D.hpp>
#include <sgpp/combigrid/pce/CombigridSurrogateModel.hpp>
#include <sgpp/combigrid/pce/CombigridSurrogateModelFactory.hpp>
#include <sgpp/combigrid/utils/AnalyticModels.hpp>
#include <sgpp/combigrid/operation/multidim/AveragingLevelManager.hpp>
#include <sgpp/combigrid/definitions.hpp>

#include <sgpp/quadrature/sampling/LatinHypercubeSampleGenerator.hpp>

#include <sgpp/globaldef.hpp>

#include <vector>

#ifdef USE_DAKOTA

BOOST_AUTO_TEST_SUITE(testPolynomialChaosExpansion)

void testPCEIshigami(std::shared_ptr<sgpp::combigrid::CombigridOperation> op,
                     std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D> basisFunction) {
  // initialize the surrogate model
  sgpp::combigrid::CombigridSurrogateModelConfiguration config;
  config.type = sgpp::combigrid::CombigridSurrogateModelsType::POLYNOMIAL_CHAOS_EXPANSION;
  config.combigridOperation = op;
  config.basisFunction = basisFunction;
  auto pce = sgpp::combigrid::createCombigridSurrogateModel(config);

  // ishigami model
  sgpp::combigrid::Ishigami ishigamiModel;

  // check the moments
  BOOST_CHECK_SMALL(std::abs(ishigamiModel.mean - pce->mean()), ishigamiModel.tolerance);
  BOOST_CHECK_SMALL(std::abs(ishigamiModel.variance - pce->variance()), ishigamiModel.tolerance);

  // check the sobol indices
  sgpp::base::DataVector sobolIndices;
  pce->getComponentSobolIndices(sobolIndices);
  for (size_t i = 0; i < sobolIndices.getSize(); i++) {
    BOOST_CHECK_SMALL(std::abs(ishigamiModel.sobolIndices[i] - sobolIndices[i]),
                      ishigamiModel.tolerance);
  }

  // check the total sobol indices
  sgpp::base::DataVector totalSobolIndices;
  pce->getTotalSobolIndices(totalSobolIndices);
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

void testPCEParbola(std::shared_ptr<sgpp::combigrid::CombigridOperation> op,
                    sgpp::combigrid::OrthogonalBasisFunctionsCollection& basisFunctions) {
  // initialize the surrogate model
  sgpp::combigrid::CombigridSurrogateModelConfiguration config;
  config.type = sgpp::combigrid::CombigridSurrogateModelsType::POLYNOMIAL_CHAOS_EXPANSION;
  config.combigridOperation = op;
  config.basisFunctions = basisFunctions;
  auto pce = sgpp::combigrid::createCombigridSurrogateModel(config);

  // check the moments
  sgpp::combigrid::Parabola parabolaModel;
  BOOST_CHECK_SMALL(std::abs(parabolaModel.mean - pce->mean()), parabolaModel.tolerance);
  BOOST_CHECK_SMALL(std::abs(parabolaModel.variance - pce->variance()), parabolaModel.tolerance);
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

  sgpp::combigrid::OrthogonalBasisFunctionsCollection functionBases;
  functionBases.push_back(functionBasis1);
  functionBases.push_back(functionBasis2);

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
    std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D> basisFunction, double mean,
    double variance, double tol) {
  // initialize the surrogate model
  sgpp::combigrid::CombigridSurrogateModelConfiguration config;
  config.type = sgpp::combigrid::CombigridSurrogateModelsType::POLYNOMIAL_STOCHASTIC_COLLOCATION;
  config.combigridOperation = op;
  config.basisFunction = basisFunction;
  auto sc = sgpp::combigrid::createCombigridSurrogateModel(config);

  // check the moments
  BOOST_CHECK_SMALL(std::abs(mean - sc->mean()), tol);
  BOOST_CHECK_SMALL(std::abs(variance - sc->variance()), tol);
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
    sgpp::combigrid::OrthogonalBasisFunctionsCollection& basisFunctions, double mean,
    double variance, double tol) {
  // initialize the surrogate model
  sgpp::combigrid::CombigridSurrogateModelConfiguration config;
  config.type = sgpp::combigrid::CombigridSurrogateModelsType::POLYNOMIAL_STOCHASTIC_COLLOCATION;
  config.combigridOperation = op;
  config.basisFunctions = basisFunctions;
  auto sc = sgpp::combigrid::createCombigridSurrogateModel(config);

  // check the moments
  BOOST_CHECK_SMALL(std::abs(mean - sc->mean()), tol);
  BOOST_CHECK_SMALL(std::abs(variance - sc->variance()), tol);
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

  sgpp::combigrid::OrthogonalBasisFunctionsCollection functionBases;
  functionBases.push_back(functionBasis1);
  functionBases.push_back(functionBasis2);

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
  auto basisFunction = std::make_shared<sgpp::combigrid::OrthogonalPolynomialBasis1D>(config);
  sgpp::combigrid::MultiFunction func(co2Model.eval);

  auto op = sgpp::combigrid::CombigridOperation::createExpClenshawCurtisPolynomialInterpolation(
      co2Model.numDims, func);
  auto op_levelManager = op->getLevelManager();

  // initialize the surrogate model
  sgpp::combigrid::CombigridSurrogateModelConfiguration surrogate_config;
  surrogate_config.type =
      sgpp::combigrid::CombigridSurrogateModelsType::POLYNOMIAL_STOCHASTIC_COLLOCATION;
  surrogate_config.combigridOperation = op;
  surrogate_config.basisFunction = basisFunction;
  auto sc = sgpp::combigrid::createCombigridSurrogateModel(surrogate_config);
  auto tensor_levelManager = sc->getConfig().combigridTensorOperation->getLevelManager();

  op_levelManager->addRegularLevels(2);
  tensor_levelManager->addLevelsFromStructure(op_levelManager->getLevelStructure());

  // check the moments
  BOOST_CHECK_SMALL(std::abs(co2Model.mean - sc->mean()), 1e-1);
  BOOST_CHECK_SMALL(std::abs(co2Model.variance - sc->variance()), 1e-1);

  auto op2 = sgpp::combigrid::CombigridOperation::createExpL2LejaPolynomialInterpolation(
      co2Model.numDims, func);
  auto op2_levelManager = op2->getLevelManager();

  // compute variance of the estimator
  sc->updateOperation(op2);
  tensor_levelManager = sc->getConfig().combigridTensorOperation->getLevelManager();

  op2_levelManager->addRegularLevels(7);
  tensor_levelManager->addLevelsFromStructure(op2_levelManager->getLevelStructure());

  // check the moments
  BOOST_CHECK_SMALL(std::abs(co2Model.mean - sc->mean()), co2Model.tolerance);
  BOOST_CHECK_SMALL(std::abs(co2Model.variance - sc->variance()), co2Model.tolerance);
}

BOOST_AUTO_TEST_CASE(testMoments) {
  // use the atan function as model function
  sgpp::combigrid::AtanBeta model;

  // -----------------------------------------------------------------------------------
  // initialize basis functions
  sgpp::combigrid::OrthogonalPolynomialBasis1DConfiguration config;
  config.polyParameters.type_ = sgpp::combigrid::OrthogonalPolynomialBasisType::JACOBI;
  config.polyParameters.lowerBound_ = model.bounds[0];
  config.polyParameters.upperBound_ = model.bounds[1];
  config.polyParameters.alpha_ = model.alpha1;
  config.polyParameters.beta_ = model.beta1;
  auto functionBasis1 = std::make_shared<sgpp::combigrid::OrthogonalPolynomialBasis1D>(config);

  config.polyParameters.type_ = sgpp::combigrid::OrthogonalPolynomialBasisType::JACOBI;
  config.polyParameters.lowerBound_ = model.bounds[0];
  config.polyParameters.upperBound_ = model.bounds[1];
  config.polyParameters.alpha_ = model.alpha2;
  config.polyParameters.beta_ = model.beta2;
  auto functionBasis2 = std::make_shared<sgpp::combigrid::OrthogonalPolynomialBasis1D>(config);

  sgpp::combigrid::OrthogonalBasisFunctionsCollection orthogFunctionBases;
  orthogFunctionBases.push_back(functionBasis1);
  orthogFunctionBases.push_back(functionBasis2);

  sgpp::combigrid::MultiFunction func(model.eval);

  // -----------------------------------------------------------------------------------
  // non-orthogonal basis function for basis transformation
  config.polyParameters.type_ = sgpp::combigrid::OrthogonalPolynomialBasisType::LEGENDRE;
  config.polyParameters.lowerBound_ = 0.0;
  config.polyParameters.upperBound_ = 1.0;
  auto legendreBasis = std::make_shared<sgpp::combigrid::OrthogonalPolynomialBasis1D>(config);

  // -----------------------------------------------------------------------------------
  // use tensor based refinement
  sgpp::combigrid::CombiEvaluators::TensorCollection tensor_evaluators(
      model.numDims, sgpp::combigrid::CombiEvaluators::tensorInterpolation(legendreBasis));
  auto tensor_op_lm = std::make_shared<sgpp::combigrid::AveragingLevelManager>();

  sgpp::combigrid::CombiHierarchies::Collection tensor_grids(
      model.numDims, sgpp::combigrid::CombiHierarchies::expLeja());
  auto tensor_op = std::make_shared<sgpp::combigrid::CombigridTensorOperation>(
      tensor_grids, tensor_evaluators, tensor_op_lm, func, true);

  tensor_op_lm->addRegularLevels(5);

  // initialize the surrogate model
  sgpp::combigrid::CombigridSurrogateModelConfiguration sc_config;
  sc_config.type = sgpp::combigrid::CombigridSurrogateModelsType::POLYNOMIAL_STOCHASTIC_COLLOCATION;
  sc_config.combigridTensorOperation = tensor_op;
  sc_config.basisFunctions = orthogFunctionBases;
  auto sc = sgpp::combigrid::createCombigridSurrogateModel(sc_config);
  // -----------------------------------------------------------------------------------
  // use tensor based refinement
  sgpp::combigrid::CombiEvaluators::TensorCollection tensor_evaluators_pce(
      model.numDims, sgpp::combigrid::CombiEvaluators::tensorInterpolation(legendreBasis));
  auto tensor_op_lm_pce = std::make_shared<sgpp::combigrid::AveragingLevelManager>();

  sgpp::combigrid::CombiHierarchies::Collection tensor_grids_pce(
      model.numDims, sgpp::combigrid::CombiHierarchies::expLeja());
  auto tensor_op_pce = std::make_shared<sgpp::combigrid::CombigridTensorOperation>(
      tensor_grids_pce, tensor_evaluators_pce, tensor_op_lm_pce, func, true);

  tensor_op_lm_pce->addLevelsFromStructure(tensor_op_lm->getLevelStructure());

  // initialize the surrogate model
  sgpp::combigrid::CombigridSurrogateModelConfiguration pce_config;
  pce_config.type = sgpp::combigrid::CombigridSurrogateModelsType::POLYNOMIAL_CHAOS_EXPANSION;
  pce_config.combigridTensorOperation = tensor_op_pce;
  pce_config.basisFunctions = orthogFunctionBases;
  auto pce = sgpp::combigrid::createCombigridSurrogateModel(pce_config);

  // -----------------------------------------------------------------------------------
  // check if the results are equal
  BOOST_CHECK_SMALL(std::abs(pce->mean() - sc->mean()), 1e-13);
  BOOST_CHECK_SMALL(std::abs(pce->variance() - sc->variance()), 1e-13);
}

BOOST_AUTO_TEST_SUITE_END()

#endif
