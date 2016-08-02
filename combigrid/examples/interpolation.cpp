// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp_combigrid.hpp>
#include <cmath>

using namespace sgpp;
using namespace sgpp::combigrid;

/**
 * The function we want to interpolate
 */
double f_2D(base::DataVector v) { return 4.0 * v[0] * v[0] * (v[1] - v[1] * v[1]); }

int main() {
  size_t numDimensions = 2;
  MultiFunction wrapper(
      f_2D);  // create a function object taking a data vector and returning a double
  auto operation =
      CombigridOperation::createLinearLejaPolynomialInterpolation(numDimensions, wrapper);

  size_t maxLevelSum = 2;
  double result = operation->evaluate(
      maxLevelSum,
      base::DataVector(std::vector<double>{
          0.5, 0.7}));  // creates levels (0, 0), (1, 0), (2, 0), (1, 1), (0, 1), (0, 2)

  // multiFunc(lambda x: x*x)
}
