// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

// include all combigrid headers
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

double f(sgpp::base::DataVector const &x) {
  double prod = 1.0;
  for (size_t i = 0; i < x.getSize(); ++i) {
    prod *= exp(-x[i]);
  }
  return prod;
}

int main() {
  size_t numDimensions = 2;
  MultiFunction wrapper(f);

  auto operation = CombigridOperation::createExpUniformBsplineInterpolation(numDimensions, wrapper);

  size_t maxLevelSum = 2;  // creates levels (0, 0), (1, 0), (2, 0), (1, 1), (0, 1), (0, 2)
  double result = operation->evaluate(maxLevelSum, DataVector(std::vector<double>{0.5, 0.7}));
  std::cout << result << " " << f(DataVector(std::vector<double>{0.5, 0.7})) << "\n";

  auto my_storage = operation->getStorage();
  std::cout << "# grid points: " << my_storage->getNumEntries() << std::endl;

  auto my_levelmanager = operation->getLevelManager();
  std::vector<DataVector> my_gridpoints = my_levelmanager->getAllGridPoints();
  for (size_t i = 0; i < my_gridpoints.size(); i++) {
    for (size_t j = 0; j < my_gridpoints[i].getSize(); j++) {
      std::cout << my_gridpoints[i][j] << " ";
    }
    std::cout << "\n";
  }
}
