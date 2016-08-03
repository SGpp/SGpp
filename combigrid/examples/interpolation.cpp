// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp_combigrid.hpp>

#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

using sgpp::base::DataVector;
using sgpp::combigrid::MultiFunction;
using sgpp::combigrid::CombigridOperation;
using sgpp::combigrid::AbstractPointHierarchy;
using sgpp::combigrid::NestedPointHierarchy;
using sgpp::combigrid::LejaPointDistribution;
using sgpp::combigrid::IdentityPointOrdering;
using sgpp::combigrid::LinearGrowthStrategy;
using sgpp::combigrid::AbstractLinearEvaluator;
using sgpp::combigrid::FloatArrayVector;
using sgpp::combigrid::FloatScalarVector;
using sgpp::combigrid::ArrayEvaluator;
using sgpp::combigrid::BarycentricInterpolationEvaluator;
using sgpp::combigrid::CombigridTreeStorage;
using sgpp::combigrid::FullGridTensorEvaluator;
using sgpp::combigrid::FunctionLookupTable;
using sgpp::combigrid::CombigridEvaluator;
using sgpp::combigrid::WeightedRatioLevelManager;

/**
 * The function we want to interpolate
 */
double f_2D(DataVector v) { return 4.0 * v[0] * v[0] * (v[1] - v[1] * v[1]); }

void simpleInterpolation() {
  size_t numDimensions = 2;
  MultiFunction wrapper(
      f_2D);  // create a function object taking a data vector and returning a double
  auto operation =
      CombigridOperation::createLinearLejaPolynomialInterpolation(numDimensions, wrapper);

  size_t maxLevelSum = 2;
  double result = operation->evaluate(
      maxLevelSum,
      DataVector(std::vector<double>{
          0.5, 0.7}));  // creates levels (0, 0), (1, 0), (2, 0), (1, 1), (0, 1), (0, 2)

  std::cout << result << "\n";
}

void adaptiveInterpolation() {
  size_t numDimensions = 2;
  MultiFunction wrapper(
      f_2D);  // create a function object taking a data vector and returning a double
  auto operation =
      CombigridOperation::createLinearLejaPolynomialInterpolation(numDimensions, wrapper);

  size_t maxLevelSum = 2;
  DataVector param(std::vector<double>{0.5, 0.7});
  auto levelManager = std::make_shared<WeightedRatioLevelManager>();

  operation->setParameters(param);
  operation->setLevelManager(levelManager);
  operation->getLevelManager()->addLevelsAdaptive(300);

  double result = operation->getResult();

  std::cout << result << "\n";
}

void multistageInterpolation() {
  size_t numDimensions = 2;

  std::vector<DataVector> gridPoints;

  auto dummyFunc = [&gridPoints](DataVector const &param) -> double {
    gridPoints.push_back(param);
    return 0.0;
  };

  std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies(numDimensions);
  size_t growthFactor = 2;  // use n = 1 + 2 * level points
  // create a point hierarchy of nested points (leja points), which are already ordered in the
  // nesting order
  pointHierarchies[0] = std::make_shared<NestedPointHierarchy>(
      std::make_shared<LejaPointDistribution>(),
      std::make_shared<IdentityPointOrdering>(std::make_shared<LinearGrowthStrategy>(growthFactor),
                                              false));  // not yet sorted
  pointHierarchies[1] = pointHierarchies[0];

  // CREATE EVALUATORS
  std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>>> evaluators(numDimensions);
  evaluators[0] = std::make_shared<ArrayEvaluator<BarycentricInterpolationEvaluator>>(
      true);  // true means that the operator needs a parameter (the interpolation point)
  evaluators[1] = evaluators[0];

  auto storage = std::make_shared<CombigridTreeStorage>(pointHierarchies, MultiFunction(dummyFunc));

  auto fullGridEval = std::make_shared<FullGridTensorEvaluator<FloatArrayVector>>(
      storage, evaluators, pointHierarchies);

  auto combiGridEval =
      std::make_shared<CombigridEvaluator<FloatArrayVector>>(numDimensions, fullGridEval);

  size_t maxLevelSum = 2;

  std::vector<FloatArrayVector> parameters(2);
  parameters[0] = FloatArrayVector(
      std::vector<FloatScalarVector>{FloatScalarVector(0.5), FloatScalarVector(0.3)});
  parameters[1] = FloatArrayVector(
      std::vector<FloatScalarVector>{FloatScalarVector(0.7), FloatScalarVector(0.27)});

  fullGridEval->setParameters(parameters);
  combiGridEval->addRegularLevels(
      maxLevelSum);  // here, the callback function inserts the grid points to gridPoints

  // ----- EVAL FUNCTION

  FunctionLookupTable funcLookup((MultiFunction(f_2D)));

  std::cout << "grid points:\n";

  for (auto &gridPoint : gridPoints) {
    std::cout << gridPoint[0] << ", " << gridPoint[1] << "\n";
    funcLookup(gridPoint);
  }

  auto realStorage = std::make_shared<CombigridTreeStorage>(
      pointHierarchies,
      MultiFunction([&funcLookup](DataVector const &param) { return funcLookup(param); }));

  auto realFullGridEval = std::make_shared<FullGridTensorEvaluator<FloatArrayVector>>(
      realStorage, evaluators, pointHierarchies);

  auto realCombiGridEval =
      std::make_shared<CombigridEvaluator<FloatArrayVector>>(numDimensions, realFullGridEval);

  realFullGridEval->setParameters(parameters);
  realCombiGridEval->addRegularLevels(maxLevelSum);

  FloatArrayVector result = realCombiGridEval->getValue();

  std::cout << "Results:\n";

  for (size_t i = 0; i < result.size(); ++i) {
    std::cout << result.get(i).getValue() << "\n";
  }
}

int main() {
  simpleInterpolation();
  adaptiveInterpolation();
  multistageInterpolation();
}
