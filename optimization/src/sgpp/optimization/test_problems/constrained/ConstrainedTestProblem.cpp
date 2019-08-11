// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/tools/RandomNumberGenerator.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/test_problems/constrained/ConstrainedTestProblem.hpp>

namespace sgpp {
namespace optimization {
namespace test_problems {

ConstrainedTestProblem::ConstrainedTestProblem(size_t d) : d(d), displacement(d, 0.0) {}

ConstrainedTestProblem::~ConstrainedTestProblem() {}

double ConstrainedTestProblem::getOptimalPoint(base::DataVector& x) {
  // reverse displace optimal point
  const double fx = getOptimalPointUndisplaced(x);
  x.sub(displacement);
  return fx;
}

void ConstrainedTestProblem::generateDisplacement() {
  generateDisplacement(DEFAULT_STANDARD_DEVIATION);
}

void ConstrainedTestProblem::generateDisplacement(double stdDev) {
  // generate displacement until a feasible one is found

  do {
    for (size_t t = 0; t < d; t++) {
      // every component is normally distributed
      displacement[t] = base::RandomNumberGenerator::getInstance().getGaussianRN(0.0, stdDev);
    }
  } while (!isDisplacementFeasible());

  // set the displacement also in the objective function and
  // in the constraint functions
  getObjectiveFunction().setDisplacement(displacement);
  getInequalityConstraintFunction().setDisplacement(displacement);
  getEqualityConstraintFunction().setDisplacement(displacement);
}

const base::DataVector& ConstrainedTestProblem::getDisplacement() const { return displacement; }

void ConstrainedTestProblem::setDisplacement(const base::DataVector& displacement) {
  this->displacement = displacement;
  // set the displacement also in the objective function and
  // in the constraint functions
  getObjectiveFunction().setDisplacement(displacement);
  getInequalityConstraintFunction().setDisplacement(displacement);
  getEqualityConstraintFunction().setDisplacement(displacement);
}

bool ConstrainedTestProblem::isDisplacementFeasible() {
  // return true, if the optimal point after displacing lies still in
  // the domain [0, 1]^d
  base::DataVector xOpt(d);
  getOptimalPoint(xOpt);

  for (size_t t = 0; t < d; t++) {
    if ((xOpt[t] < 0.0) || (xOpt[t] > 1.0)) {
      return false;
    }
  }

  return true;
}
}  // namespace test_problems
}  // namespace optimization
}  // namespace sgpp
