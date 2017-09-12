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
  //  size_t numDimensions = 2;
  size_t numDimensions = 1;
  //  MultiFunction wrapper(f_2D);
  MultiFunction wrapper(f_1D);
  auto operation = CombigridOperation::createExpUniformBsplineInterpolation(numDimensions, wrapper);

  size_t maxLevelSum = 3;
  //  double result = operation->evaluate(maxLevelSum, DataVector(std::vector<double>{0.5, 0.5}));
  double result = operation->evaluate(maxLevelSum, DataVector(std::vector<double>{0.771}));

  std::cout << "result: " << result << " f(x)" << f_1D(DataVector(std::vector<double>{0.771}))
            << "\n";
  //
  //  auto levelmanager = operation->getLevelManager();
  //  std::vector<sgpp::base::DataVector> gridpoints = levelmanager->getAllGridPoints();

  // test evaluating a Bspline
  //  size_t degree = 3;
  //  std::cout << "------------\n\n";
  //  std::unique_ptr<sgpp::base::SBsplineBase> bsplineBasis;
  //  //  bsplineBasis = std::unique_ptr<sgpp::base::SBsplineBase>(new
  //  //  sgpp::base::SBsplineBase(degree));
  //  double level = 4;
  //  double k = 1;
  //  std::vector<double> points = {0, 0.0625, 0.125, 0.1875, 0.25, 0.3125, 0.375, 0.4375, 0.5};
  //  for (unsigned int i = 0; i < points.size(); i++) {
  //    //    std::cout << bsplineBasis->eval(1, 1, points[i]) << " ";
  //    std::cout << bsplineBasis->uniformBSpline(points[i] * pow(2, level) - k, degree) << " ";
  //}
}
