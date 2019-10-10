// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/test_problems/TestVectorFunction.hpp>

namespace sgpp {
namespace optimization {
namespace test_problems {

TestVectorFunction::TestVectorFunction(size_t d, size_t m)
    : base::VectorFunction(d, m), displacement(d, 0.0), xTmp(d) {}

TestVectorFunction::~TestVectorFunction() {}

void TestVectorFunction::eval(const base::DataVector& x, base::DataVector& value) {
  // displace vector before evaluation
  for (size_t t = 0; t < d; t++) {
    xTmp[t] = x[t] + displacement[t];
  }

  evalUndisplaced(xTmp, value);
}

const base::DataVector& TestVectorFunction::getDisplacement() const { return displacement; }

void TestVectorFunction::setDisplacement(const base::DataVector& displacement) {
  this->displacement = displacement;
}
}  // namespace test_problems
}  // namespace optimization
}  // namespace sgpp
