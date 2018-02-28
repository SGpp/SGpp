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
#include <sgpp/combigrid/algebraic/FirstMomentNormStrategy.hpp>
#include <sgpp/combigrid/algebraic/VarianceNormStrategy.hpp>
#include <sgpp/combigrid/GeneralFunction.hpp>

#include <sgpp/combigrid/definitions.hpp>
#include <sgpp/globaldef.hpp>

#include <vector>

#ifdef USE_DAKOTA

BOOST_AUTO_TEST_SUITE(testRefinement)

BOOST_AUTO_TEST_CASE(testVarianceBasedRefinement) {
  // use the atan function as model function
  sgpp::combigrid::AtanUniform model;
  sgpp::combigrid::MultiFunction func(model.eval);

  size_t regularLevel = 1;

  // -----------------------------------------------------------------------------------
  // use tensor based refinement
  sgpp::combigrid::OrthogonalPolynomialBasis1DConfiguration config;
  config.polyParameters.type_ = sgpp::combigrid::OrthogonalPolynomialBasisType::LEGENDRE;
  auto functionBasis = std::make_shared<sgpp::combigrid::OrthogonalPolynomialBasis1D>(config);

  sgpp::combigrid::CombiEvaluators::TensorCollection tensor_evaluators(
      model.numDims, sgpp::combigrid::CombiEvaluators::tensorInterpolation(functionBasis));
  auto tensor_op_lm = std::make_shared<sgpp::combigrid::AveragingLevelManager>();

  sgpp::combigrid::CombiHierarchies::Collection tensor_grids(
      model.numDims, sgpp::combigrid::CombiHierarchies::expLeja());
  auto tensor_op = std::make_shared<sgpp::combigrid::CombigridTensorOperation>(
      tensor_grids, tensor_evaluators, tensor_op_lm, func, true,
      sgpp::combigrid::FullGridSummationStrategyType::TENSORVARIANCE);

  //  std::cout << "------------------------------------------------------" << std::endl;
  //  std::cout << "Tensor" << std::endl;
  tensor_op_lm->addRegularLevels(regularLevel);
  tensor_op_lm->addLevelsAdaptive(50);
  //  tensor_op_lm->addLevelsAdaptive(100);
  //  tensor_op_lm->addLevelsAdaptive(150);

  // -----------------------------------------------------------------------------------
  // use quadrature based refinement
  sgpp::combigrid::CombiEvaluators::MultiCollection evaluators(
      model.numDims, sgpp::combigrid::CombiEvaluators::polynomialScalarProduct());
  auto variance_op_lm = std::make_shared<sgpp::combigrid::AveragingLevelManager>();

  sgpp::combigrid::CombiHierarchies::Collection grids(model.numDims,
                                                      sgpp::combigrid::CombiHierarchies::expLeja());
  auto variance_op = std::make_shared<sgpp::combigrid::CombigridMultiOperation>(
      grids, evaluators, variance_op_lm, func, true,
      sgpp::combigrid::FullGridSummationStrategyType::VARIANCE);

  //  std::cout << "------------------------------------------------------" << std::endl;
  //  std::cout << "Quadrature" << std::endl;
  variance_op_lm->addRegularLevels(regularLevel);
  variance_op_lm->addLevelsAdaptive(50);
  //  variance_op_lm->addLevelsAdaptive(100);
  //  variance_op_lm->addLevelsAdaptive(150);

  // check if both approaches lead to the same result
  //  std::cout << "num grid points: " << tensor_op->numGridPoints() << " = "
  //            << variance_op_lm->numGridPoints() << std::endl;
  BOOST_CHECK_EQUAL(tensor_op->numGridPoints(), variance_op_lm->numGridPoints());

  // check the whole level structure
  auto tensor_infos = tensor_op_lm->getInfoOnAddedLevels()->getInfos();
  auto quadrature_infos = variance_op_lm->getInfoOnAddedLevels()->getInfos();

  //  std::cout << "num refinement steps: " << tensor_infos->size() << " = " <<
  //  quadrature_infos->size()
  //            << std::endl;
  BOOST_CHECK_EQUAL(tensor_infos->size(), quadrature_infos->size());

  for (size_t i = 0; i < tensor_infos->size(); i++) {
    auto tensor_info_map = (*tensor_infos)[i];
    auto quadrature_info_map = (*quadrature_infos)[i];

    //    std::cout << "#i=" << i << ": num new subspaces: " << tensor_info_map->size() << " = "
    //              << quadrature_info_map->size() << std::endl;
    BOOST_CHECK_EQUAL(tensor_info_map.size(), quadrature_info_map.size());

    for (auto tensor_info_iterator : tensor_info_map) {
      auto level = tensor_info_iterator.first;

      //      std::cout << "  (" << level[0] << ", " << level[1] << ") -> "
      //                << (quadrature_info_map->find(level) != quadrature_info_map->end());
      BOOST_CHECK(quadrature_info_map.find(level) != quadrature_info_map.end());
      auto tensor_level_info = tensor_info_iterator.second;
      auto quadrature_level_info = quadrature_info_map.find(level)->second;

      // check if the same difference norms where computed
      //      std::cout << "; priority: " << tensor_level_info->priority << " = "
      //                << quadrature_level_info->priority << "; norms: " << tensor_level_info->norm
      //                << " = " << quadrature_level_info->norm << std::endl;
      BOOST_CHECK_SMALL(std::abs(tensor_level_info.priority - quadrature_level_info.priority),
                        1e-13);
      BOOST_CHECK_SMALL(std::abs(tensor_level_info.norm - quadrature_level_info.norm), 1e-13);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()

#endif
