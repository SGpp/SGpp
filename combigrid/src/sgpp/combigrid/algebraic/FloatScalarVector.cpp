// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/algebraic/FloatScalarVector.hpp>

#include <iostream>

namespace sgpp {
namespace combigrid {
std::ostream &operator<<(std::ostream &stream, FloatScalarVector v) {
  stream << v.getValue();
  return stream;
}

std::istream &operator>>(std::istream &stream, FloatScalarVector &v) {
  double val = 0.0;
  stream >> val;
  v = FloatScalarVector(val);
  return stream;
}

} /* namespace combigrid */
} /* namespace sgpp*/
