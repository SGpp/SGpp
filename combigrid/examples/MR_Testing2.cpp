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

double f_1D(DataVector v) { return 4.0 * v[0] * v[0]; }

double f_2D(DataVector v) { return 4.0 * v[0] * v[0] * (v[1] - v[1] * v[1]); }

typedef sgpp::combigrid::AveragingLevelManager StandardLevelManager;

int main() {
  //  size_t numDimensions = 2;
  size_t numDimensions = 1;
  //  MultiFunction wrapper(f_2D);
  MultiFunction wrapper(f_1D);
  auto operation = CombigridOperation::createExpUniformBsplineInterpolation(numDimensions, wrapper);

  size_t maxLevelSum = 3;
  //  double result = operation->evaluate(maxLevelSum, DataVector(std::vector<double>{0.5, 0.5}));
  double result = operation->evaluate(maxLevelSum, DataVector(std::vector<double>{0.5}));
  result = operation->evaluate(maxLevelSum, DataVector(std::vector<double>{0.25}));
  result = operation->evaluate(maxLevelSum, DataVector(std::vector<double>{0.75}));

  std::cout << result << "\n";
  //
  //  auto levelmanager = operation->getLevelManager();
  //  std::vector<sgpp::base::DataVector> gridpoints = levelmanager->getAllGridPoints();
}
