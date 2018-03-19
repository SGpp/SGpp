// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <sgpp/combigrid/grid/distribution/ClenshawCurtisDistribution.hpp>
#include <sgpp/combigrid/grid/distribution/LejaPointDistribution.hpp>
#include <sgpp/combigrid/grid/distribution/UniformBoundaryPointDistribution.hpp>
#include <sgpp/combigrid/grid/growth/ExponentialGrowthStrategy.hpp>
#include <sgpp/combigrid/grid/growth/LinearGrowthStrategy.hpp>
#include <sgpp/combigrid/grid/hierarchy/NestedPointHierarchy.hpp>
#include <sgpp/combigrid/grid/hierarchy/NonNestedPointHierarchy.hpp>
#include <sgpp/combigrid/grid/ordering/ExponentialLevelorderPointOrdering.hpp>
#include <sgpp/combigrid/grid/ordering/IdentityPointOrdering.hpp>
#include <sgpp/combigrid/integration/MCIntegrator.hpp>
#include <sgpp/combigrid/operation/CombigridMultiOperation.hpp>
#include <sgpp/combigrid/operation/CombigridOperation.hpp>
#include <sgpp/combigrid/operation/CombigridTensorOperation.hpp>
#include <sgpp/combigrid/operation/multidim/AveragingLevelManager.hpp>
#include <sgpp/combigrid/operation/multidim/RegularLevelManager.hpp>
#include <sgpp/combigrid/operation/multidim/CombigridEvaluator.hpp>
#include <sgpp/combigrid/operation/multidim/fullgrid/FullGridCallbackEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/ArrayEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/LinearInterpolationEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/PolynomialInterpolationEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/PolynomialQuadratureEvaluator.hpp>
#include <sgpp/combigrid/functions/OrthogonalPolynomialBasis1D.hpp>
#include <sgpp/combigrid/storage/tree/CombigridTreeStorage.hpp>
#include <sgpp/combigrid/utils/Stopwatch.hpp>
#include <sgpp/combigrid/utils/Utils.hpp>
#include <sgpp/combigrid/utils/AnalyticModels.hpp>

#include <sgpp/globaldef.hpp>

#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>
#include <iomanip>

using sgpp::base::DataVector;
using sgpp::combigrid::AbstractMultiStorage;
using sgpp::combigrid::FloatArrayVector;
using sgpp::combigrid::CombigridMultiOperation;
using sgpp::combigrid::MultiFunction;
using sgpp::combigrid::Stopwatch;
using sgpp::combigrid::MCIntegrator;
using sgpp::combigrid::AbstractPointHierarchy;
using sgpp::combigrid::CombigridEvaluator;
using sgpp::combigrid::FloatScalarVector;
using sgpp::combigrid::NestedPointHierarchy;
using sgpp::combigrid::LejaPointDistribution;
using sgpp::combigrid::IdentityPointOrdering;
using sgpp::combigrid::LinearGrowthStrategy;
using sgpp::combigrid::AbstractLinearEvaluator;
using sgpp::combigrid::PolynomialInterpolationEvaluator;
using sgpp::combigrid::CombigridTreeStorage;
using sgpp::combigrid::AveragingLevelManager;

double testFunction(DataVector const &coordinates);

double testFunction2(DataVector const &coordinates);

double testFunction3(DataVector const &coordinates);

double testFunction4(DataVector const &coordinates);

double testFunction5(DataVector const &x);

double testFunction6(DataVector const &x);

double testFunction7(DataVector const &x);

double testFunctionAtan(DataVector const &x);

/*
 void printCTResults(size_t d, size_t q) {
 const size_t samples = 10;
 auto ctInterpolator = CombigridOperation::createExpClenshawCurtisPolynomialInterpolation(d,
 testFunction);
 auto domain = std::vector<std::pair<double, double>>(d, std::pair<double, double>(0.0,
 1.0));

 MCIntegrator integrator([&](DataVector const &x) -> double {
 double diff = testFunction(x) - ctInterpolator->evaluate(q, x);
 return diff * diff;
 });

 std::cout << "d = " << d << ", q = " << q << ": " << std::sqrt(integrator.average(domain,
 samples)) << std::endl;
 }
 */

BOOST_AUTO_TEST_CASE(testLevelManagerParallel) {
  size_t numDimensions = 2;

  auto func = testFunction2;

  size_t growthFactor = 2;  // use n = 1 + 2 * level points
  std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies(
      numDimensions, std::make_shared<NestedPointHierarchy>(
                         std::make_shared<LejaPointDistribution>(),
                         std::make_shared<IdentityPointOrdering>(
                             std::make_shared<LinearGrowthStrategy>(growthFactor), false)));
  // create a point hierarchy of nested points (leja points), which are already ordered in the
  // nesting order

  // CREATE EVALUATORS
  std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>> evaluators(
      numDimensions, std::make_shared<PolynomialInterpolationEvaluator>());

  auto storage = std::make_shared<CombigridTreeStorage>(pointHierarchies, MultiFunction(func));

  auto fullGridEval =
      std::make_shared<sgpp::combigrid::FullGridCallbackEvaluator<FloatScalarVector>>(
          storage, evaluators, pointHierarchies);

  auto combiGridEval =
      std::make_shared<CombigridEvaluator<FloatScalarVector>>(numDimensions, fullGridEval);

  auto levelManager = std::make_shared<AveragingLevelManager>(combiGridEval);

  size_t maxLevelSum = 3;

  std::vector<FloatScalarVector> parameters(2);
  parameters[0] = FloatScalarVector(0.378934);
  parameters[1] = FloatScalarVector(0.89340273);

  fullGridEval->setParameters(parameters);
  levelManager->addRegularLevelsParallel(maxLevelSum, 4);
  levelManager->addLevelsAdaptiveParallel(300, 4);

  // std::cout << "test_level_manager: ";
  // std::cout << std::abs(combiGridEval->getValue().getValue() -
  //                       func(DataVector(std::vector<double>{0.378934, 0.89340273})))
  //           << "\n";
}

BOOST_AUTO_TEST_CASE(testLevelManagerAdaptive) {
  size_t numDims = 6;
  sgpp::combigrid::Genz model;
  sgpp::combigrid::MultiFunction func(model.eval);
  auto op = sgpp::combigrid::CombigridOperation::createExpClenshawCurtisPolynomialInterpolation(
      numDims, func);

  auto levelManager = std::make_shared<AveragingLevelManager>();
  op->setLevelManager(levelManager);
  levelManager->disableStatsCollection();
  size_t i = 0;
  size_t maxIterations = 30;
  while (i < maxIterations) {
    levelManager->addLevelsAdaptiveByNumLevels(10);
    BOOST_CHECK_EQUAL(op->getUpperPointBound(), op->numGridPoints());
    i++;
  }
}

#ifdef USE_DAKOTA

BOOST_AUTO_TEST_CASE(testLevelManagerStats) {
  auto func = MultiFunction(testFunctionAtan);
  size_t d = 2;
  sgpp::combigrid::OrthogonalPolynomialBasis1DConfiguration config;
  config.polyParameters.type_ = sgpp::combigrid::OrthogonalPolynomialBasisType::LEGENDRE;
  auto functionBasis = std::make_shared<sgpp::combigrid::OrthogonalPolynomialBasis1D>(config);

  auto op = sgpp::combigrid::CombigridTensorOperation::createExpLejaPolynomialInterpolation(
      functionBasis, d, func);
  auto levelManager = op->getLevelManager();
  levelManager->enableStatsCollection();

  // add regular levels
  size_t maxLevel = 8;
  for (size_t level = 0; level < maxLevel; level++) {
    levelManager->addRegularLevels(level);
  }

  // reference stats
  std::vector<double> refStats{4.00904112129312,   1.23827837473378, 0.771152084553602,
                               0.321442518162072,  0.2398746901928,  0.143988691577309,
                               0.0775274779950517, 0.013446805563535};

  auto stats = levelManager->getInfoOnAddedLevels();

  // evaluate stats
  // size_t i = 0;
  // for (auto &istats : *stats->getInfos()) {
  //   std::cout << " - - - - - - - - - - - - " << std::endl;
  //   std::cout << "iteration = " << ++i << std::endl;
  //   for (auto &item : istats) {
  //     auto &level = item.first;
  //     auto &levelInfo = item.second;
  //     std::cout << level[0] << ", " << level[1] << " = " << std::setprecision(15) <<
  //     levelInfo.norm
  //               << ", " << levelInfo.priority << std::endl;
  //   }
  // }

  // compute maximum norm per iteration
  sgpp::base::DataVector maxNorms;
  stats->maxNormPerIteration(maxNorms);
  for (size_t i = 0; i < maxNorms.getSize(); i++) {
    BOOST_CHECK_SMALL(std::abs(refStats[i] - maxNorms[i]), 1e-5);
  }
}

BOOST_AUTO_TEST_CASE(testLevelManagerStatsConversion) {
  auto func = MultiFunction(testFunctionAtan);
  size_t d = 2;

  auto op = sgpp::combigrid::CombigridOperation::createExpLejaPolynomialInterpolation(d, func);
  op->getLevelManager()->enableStatsCollection();

  sgpp::combigrid::OrthogonalPolynomialBasis1DConfiguration config;
  config.polyParameters.type_ = sgpp::combigrid::OrthogonalPolynomialBasisType::LEGENDRE;
  auto functionBasis = std::make_shared<sgpp::combigrid::OrthogonalPolynomialBasis1D>(config);

  // copy level structure to tensor grid
  auto tensor_op = sgpp::combigrid::CombigridTensorOperation::createExpLejaPolynomialInterpolation(
      functionBasis, d, func);
  tensor_op->getLevelManager()->enableStatsCollection();

  // add regular levels
  size_t maxLevel = 8;
  for (size_t level = 0; level < maxLevel; level++) {
    op->getLevelManager()->addRegularLevels(level);
    tensor_op->getLevelManager()->addLevelsFromStructure(
        op->getLevelManager()->getLevelStructure());
  }

  // evaluate stats
  //  size_t i = 0;
  //  for (auto &istats : *op->getLevelManager()->getInfoOnAddedLevels()->getInfos()) {
  //    std::cout << " - - - - - - - - - - - - " << std::endl;
  //    std::cout << "iteration = " << ++i << std::endl;
  //    for (auto &item : *istats) {
  //      auto &level = item.first;
  //      auto &levelInfo = item.second;
  //      std::cout << level[0] << ", " << level[1] << " = " << std::setprecision(15) <<
  //      levelInfo->norm
  //                << ", " << levelInfo->priority << std::endl;
  //    }
  //  }

  // reference stats
  std::vector<double> refStats{4.00904112129312,   1.23827837473378, 0.771152084553602,
                               0.321442518162072,  0.2398746901928,  0.143988691577309,
                               0.0775274779950517, 0.013446805563535};

  // check if the stats of the converted operation are correct
  sgpp::base::DataVector maxNorms;
  tensor_op->getLevelManager()->getInfoOnAddedLevels()->maxNormPerIteration(maxNorms);
  for (size_t i = 0; i < maxNorms.getSize(); i++) {
    BOOST_CHECK_SMALL(std::abs(refStats[i] - maxNorms[i]), 1e-5);
  }
}

#endif
