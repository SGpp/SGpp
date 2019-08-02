// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/test_problems/TestScalarFunction.hpp>

namespace sgpp {
namespace optimization {
namespace test_problems {

TestScalarFunction::TestScalarFunction(size_t d)
    : base::ScalarFunction(d), displacement(d, 0.0), xTmp(d) {}

TestScalarFunction::~TestScalarFunction() {}

double TestScalarFunction::eval(const base::DataVector& x) {
  // displace vector before evaluation
  for (size_t t = 0; t < d; t++) {
    xTmp[t] = x[t] + displacement[t];
  }

  return evalUndisplaced(xTmp);
}

const base::DataVector& TestScalarFunction::getDisplacement() const { return displacement; }

void TestScalarFunction::setDisplacement(const base::DataVector& displacement) {
  this->displacement = displacement;
}
}  // namespace test_problems
}  // namespace optimization
}  // namespace sgpp
