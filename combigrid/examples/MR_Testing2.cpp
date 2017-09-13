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

#include <sgpp/base/operation/hash/common/basis/BsplineBasis.hpp>

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
  //  size_t numDimensions = 1;
  size_t numDimensions = 2;
  //  MultiFunction wrapper(f_1D);
  MultiFunction wrapper(f_2D);

  auto operation = CombigridOperation::createExpUniformBsplineInterpolation(numDimensions, wrapper);
  //  auto operation = CombigridOperation::createExpUniformLinearInterpolation(numDimensions,
  //  wrapper);

  size_t maxLevelSum = 2;
  //  double result = operation->evaluate(maxLevelSum, DataVector(std::vector<double>{0.771}));
  double result = operation->evaluate(maxLevelSum, DataVector(std::vector<double>{0.5, 0.5}));

  std::cout << "result: " << result << " f(x)=" << f_2D(DataVector(std::vector<double>{0.5, 0.5}))
            << "\n";
  //
  //  auto levelmanager = operation->getLevelManager();
  //  std::vector<sgpp::base::DataVector> gridpoints = levelmanager->getAllGridPoints();
}
