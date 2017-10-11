// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/combigrid/functions/AbstractInfiniteFunctionBasis1D.hpp>

namespace sgpp {
namespace combigrid {

class MonomialFunctionBasis1D : public AbstractInfiniteFunctionBasis1D {
 public:
  ~MonomialFunctionBasis1D();

  virtual double evaluate(size_t basisIndex, double xValue);
};

} /* namespace combigrid */
} /* namespace sgpp */
