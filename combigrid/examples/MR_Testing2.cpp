// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

// include all combigrid headers
#include <sgpp/combigrid/operation/Configurations.hpp>
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
using sgpp::combigrid::PolynomialInterpolationEvaluator;
using sgpp::combigrid::CombigridTreeStorage;
using sgpp::combigrid::FullGridLinearCallbackEvaluator;
using sgpp::combigrid::FunctionLookupTable;
using sgpp::combigrid::CombigridEvaluator;
using sgpp::combigrid::WeightedRatioLevelManager;

double f_2D(DataVector v) { return 4.0 * v[0] * v[0] * (v[1] - v[1] * v[1]); }

typedef sgpp::combigrid::AveragingLevelManager StandardLevelManager;

void testInterpolationWithBsplines() {
  size_t numDimensions = 1;
  sgpp::combigrid::LinearInterpolationEvaluator eval;
  auto evaluator = std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>>(
      numDimensions, sgpp::combigrid::CombiEvaluators::linearInterpolation());

  auto pointhierarchy = std::vector<std::shared_ptr<AbstractPointHierarchy>>(
      numDimensions, sgpp::combigrid::CombiHierarchies::expUniform());

  auto levelmanager = std::make_shared<StandardLevelManager>();

  std::vector<double> points;
  points.push_back(0.0);
  points.push_back(0.5);
  points.push_back(1.0);

  eval.setGridPoints(points);
}

int main() {
  size_t numDimensions = 2;
  MultiFunction wrapper(f_2D);
  auto operation = CombigridOperation::createExpUniformBsplineInterpolation(numDimensions, wrapper);

  size_t maxLevelSum = 2;
  double result = operation->evaluate(maxLevelSum, DataVector(std::vector<double>{0.5, 0.7}));

  std::cout << result << "\n";
  //
  //  auto levelmanager = operation->getLevelManager();
  //  std::vector<sgpp::base::DataVector> gridpoints = levelmanager->getAllGridPoints();
}
