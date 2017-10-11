// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/combigrid/functions/AbstractInfiniteFunctionBasis1D.hpp>

#ifdef USE_DAKOTA
#include <BasisPolynomial.hpp>

namespace sgpp {
namespace combigrid {

class LegendreBasis1D : public AbstractInfiniteFunctionBasis1D {
 public:
  LegendreBasis1D();
  virtual ~LegendreBasis1D();

  double evaluate(size_t basisIndex, double xValue) override;

  Pecos::BasisPolynomial basisPoly;
};

} /* namespace combigrid */
} /* namespace sgpp */

#endif
