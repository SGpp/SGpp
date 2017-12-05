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
#include <sgpp/combigrid/operation/multidim/AveragingLevelManager.hpp>
#include <sgpp/combigrid/operation/Configurations.hpp>

#include <sgpp/combigrid/definitions.hpp>
#include <sgpp/globaldef.hpp>

#include <vector>

#ifdef USE_DAKOTA

BOOST_AUTO_TEST_SUITE(testRefinement)

BOOST_AUTO_TEST_CASE(testVarianceBasedRefinement) {
  // use the ishigami function as model function
  sgpp::combigrid::Ishigami ishigamiModel;
  sgpp::combigrid::MultiFunction func(ishigamiModel.eval);

  size_t regularLevel = 1;
  size_t numAdaptivePoints = 200;

  // use tensor based refinement
  sgpp::combigrid::OrthogonalPolynomialBasis1DConfiguration config;
  config.polyParameters.type_ = sgpp::combigrid::OrthogonalPolynomialBasisType::LEGENDRE;
  auto functionBasis = std::make_shared<sgpp::combigrid::OrthogonalPolynomialBasis1D>(config);

  auto tensor_op =
      sgpp::combigrid::CombigridTensorOperation::createExpL2LejaPolynomialInterpolation(
          functionBasis, ishigamiModel.numDims, func);
  auto tensor_op_lm = tensor_op->getLevelManager();
  tensor_op_lm->addRegularLevels(regularLevel);
  tensor_op_lm->addLevelsAdaptive(numAdaptivePoints);

  // use quadrature based refinement
  sgpp::combigrid::CombiEvaluators::MultiCollection evaluators(
      ishigamiModel.numDims, sgpp::combigrid::CombiEvaluators::polynomialScalarProduct());
  auto levelManager = std::make_shared<sgpp::combigrid::AveragingLevelManager>();

  sgpp::combigrid::CombiHierarchies::Collection grids(
      ishigamiModel.numDims, sgpp::combigrid::CombiHierarchies::expL2Leja());
  auto variance_op = std::make_shared<sgpp::combigrid::CombigridMultiOperation>(
      grids, evaluators, levelManager, func, true,
      sgpp::combigrid::FullGridSummationStrategyType::VARIANCE);

  auto variance_op_lm = variance_op->getLevelManager();
  variance_op_lm->addRegularLevels(regularLevel);
  variance_op_lm->addLevelsAdaptive(numAdaptivePoints);

  // check if both approaches lead to the same result
  std::cout << tensor_op->numGridPoints() << " = " << variance_op_lm->numGridPoints() << std::endl;
}

BOOST_AUTO_TEST_SUITE_END()

#endif
