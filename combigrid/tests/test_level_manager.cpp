// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <sgpp/combigrid/grid/distribution/ClenshawCurtisDistribution.hpp>
#include <sgpp/combigrid/grid/distribution/LejaPointDistribution.hpp>
#include <sgpp/combigrid/grid/distribution/UniformPointDistribution.hpp>
#include <sgpp/combigrid/grid/growth/ExponentialGrowthStrategy.hpp>
#include <sgpp/combigrid/grid/growth/LinearGrowthStrategy.hpp>
#include <sgpp/combigrid/grid/hierarchy/NestedPointHierarchy.hpp>
#include <sgpp/combigrid/grid/hierarchy/NonNestedPointHierarchy.hpp>
#include <sgpp/combigrid/grid/ordering/ExponentialLevelorderPointOrdering.hpp>
#include <sgpp/combigrid/grid/ordering/IdentityPointOrdering.hpp>
#include <sgpp/combigrid/integration/MCIntegrator.hpp>
#include <sgpp/combigrid/operation/CombigridMultiOperation.hpp>
#include <sgpp/combigrid/operation/CombigridOperation.hpp>
#include <sgpp/combigrid/operation/multidim/AveragingLevelManager.hpp>
#include <sgpp/combigrid/operation/multidim/CombigridEvaluator.hpp>
#include <sgpp/combigrid/operation/multidim/FullGridTensorEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/ArrayEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/BarycentricInterpolationEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/LinearInterpolationEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/QuadratureEvaluator.hpp>
#include <sgpp/combigrid/storage/tree/CombigridTreeStorage.hpp>
#include <sgpp/combigrid/utils/Stopwatch.hpp>
#include <sgpp/combigrid/utils/Utils.hpp>
#include <sgpp/globaldef.hpp>

#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

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
using sgpp::combigrid::BarycentricInterpolationEvaluator;
using sgpp::combigrid::FullGridTensorEvaluator;
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
      numDimensions, std::make_shared<BarycentricInterpolationEvaluator>());

  auto storage = std::make_shared<CombigridTreeStorage>(pointHierarchies, MultiFunction(func));

  auto fullGridEval = std::make_shared<FullGridTensorEvaluator<FloatScalarVector>>(
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

  std::cout << "test_level_manager: ";
  std::cout << std::abs(combiGridEval->getValue().getValue() -
                        func(DataVector(std::vector<double>{0.378934, 0.89340273})))
            << "\n";
}
