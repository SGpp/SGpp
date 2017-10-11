// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/functions/LegendreBasis1D.hpp>

#ifdef USE_DAKOTA
#include <BasisPolynomial.hpp>

namespace sgpp {
namespace combigrid {

LegendreBasis1D::LegendreBasis1D() : basisPoly(Pecos::LEGENDRE_ORTHOG) {}

LegendreBasis1D::~LegendreBasis1D() {}

double LegendreBasis1D::evaluate(size_t basisIndex, double xValue) {
  double invNorm = 1. / std::sqrt(basisPoly.norm_squared(basisIndex));
  return invNorm * basisPoly.type1_value(xValue, basisIndex);
}

} /* namespace combigrid */
} /* namespace sgpp */

#endif
