// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

/**
 * First, some includes.
 */

#include <sgpp/combigrid/operation/CombigridMultiOperation.hpp>
#include <sgpp/combigrid/operation/CombigridOperation.hpp>

#include <cmath>

#include <iostream>

/**
 * The first thing we need is a function to evaluate. This function will be evaluated on the domain
 * [0, 1]^d. This particular function can be used with any number of dimensions.
 */

double f(sgpp::base::DataVector const &x) {
  double prod = 1.0;
  for (size_t i = 0; i < x.getSize(); ++i) {
    prod *= exp(-x[i]);
  }
  return prod;
}

// We have to wrap f in a sgpp::combigrid::MultiFunction object.
static const sgpp::combigrid::MultiFunction func(f);

// Let's use a 3D-function.
size_t d = 3;

/**
 * Here comes the first and very simple example.
 */
void example1() {
  // Let's increase the number of points by two for each level.
  size_t growthFactor = 2;

  // Now create the operation object that handles the evaluation. The evaluation mode is quadrature,
  // so it will approximate the integral of f over [0, 1]^d. It uses Leja points with 1 + 2*l
  // points in level l. The level starts from zero, higher level means finer grid.
  // Slower growth of the number of points per level means that the total number of points used can
  // be controlled better.
  std::shared_ptr<sgpp::combigrid::CombigridOperation> operation =
      sgpp::combigrid::CombigridOperation::createLinearLejaQuadrature(d, func, growthFactor);

  // Now, we compute the result. The parameter 2 means that grid at level-multi-indices with a
  // 1-norm (i.e. sum of entries) less than or equal to 2 are used. In our 3D case, these are
  // exactly the levels (0, 0, 0), (1, 0, 0), (2, 0, 0), (1, 1, 0), (1, 0, 1), (0, 1, 0), (0, 2, 0),
  // (0, 1, 1), (0, 0, 1) and (0, 0, 2).
  double result = operation->evaluate(2);

  // Now compare the result to the analytical solution:
  std::cout << "Quadrature result: " << result
            << ", analytical solution: " << pow(1 - 1.0 / M_E, static_cast<double>(d)) << "\n";

  // We can also find out how many function evaluations have been used by accessing the storage
  // which stores computed function values:
  std::cout << "Number of function evaluations: " << operation->getStorage()->getNumEntries()
            << "\n";
}

void example2() {}

int main() {
  example1();
  example2();
}
