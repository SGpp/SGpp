// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <cstddef>

namespace sgpp {
namespace combigrid {

class AbstractInfiniteFunctionBasis1D {
 public:
  AbstractInfiniteFunctionBasis1D();
  virtual ~AbstractInfiniteFunctionBasis1D();

  virtual double evaluate(size_t basisIndex, double xValue) = 0;
};

} /* namespace combigrid */
} /* namespace sgpp */
