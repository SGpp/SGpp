// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/functions/LegendreBasis1D.hpp>

#include <iostream>

#ifdef USE_DAKOTA
#include <BasisPolynomial.hpp>
#endif

namespace sgpp {
namespace combigrid {

#ifdef USE_DAKOTA
LegendreBasis1D::LegendreBasis1D() : basisPoly(Pecos::LEGENDRE_ORTHOG) {}
#else
LegendreBasis1D::LegendreBasis1D() {}
#endif

LegendreBasis1D::~LegendreBasis1D() {}

double LegendreBasis1D::evaluate(size_t basisIndex, double xValue) {
#ifdef USE_DAKOTA
  double invNorm = 1. / std::sqrt(basisPoly.norm_squared(static_cast<short unsigned>(basisIndex)));
  return invNorm *
         basisPoly.type1_value(2.0 * xValue - 1.0, static_cast<short unsigned>(basisIndex));
#else
  std::cerr << "Error in LegendreBasis1D::evaluate: "
            << "SG++ was compiled without DAKOTA support!\n";
  return false;
#endif
}

} /* namespace combigrid */
} /* namespace sgpp */
