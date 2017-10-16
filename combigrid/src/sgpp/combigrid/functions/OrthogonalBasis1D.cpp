// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/functions/OrthogonalBasis1D.hpp>

#include <iostream>

#ifdef USE_DAKOTA
#include <BasisPolynomial.hpp>
#include <pecos_global_defs.hpp>
#endif

#ifndef M_SQRT2
#define M_SQRT2 1.41421356237309504880
#endif

namespace sgpp {
namespace combigrid {

OrthogonalBasis1D::OrthogonalBasis1D(OrthogonalPolynomialBasisType basisType)
    : basisType(basisType) {
#ifdef USE_DAKOTA
  switch (basisType) {
    case OrthogonalPolynomialBasisType::LEGENDRE:
      basisPoly = std::make_shared<Pecos::BasisPolynomial>(Pecos::LEGENDRE_ORTHOG);
      break;
    case OrthogonalPolynomialBasisType::HERMITE:
      basisPoly = std::make_shared<Pecos::BasisPolynomial>(Pecos::HERMITE_ORTHOG);
      break;
    case OrthogonalPolynomialBasisType::JACOBI:
      basisPoly = std::make_shared<Pecos::BasisPolynomial>(Pecos::JACOBI_ORTHOG);
      break;
  }
#endif
}

OrthogonalBasis1D::~OrthogonalBasis1D() {}

double OrthogonalBasis1D::normalizeInput(double xValue) {
  switch (basisType) {
    case OrthogonalPolynomialBasisType::LEGENDRE:
    case OrthogonalPolynomialBasisType::JACOBI:
      // [0, 1] -> [-1, 1]
      return 2.0 * xValue - 1.0;
    case OrthogonalPolynomialBasisType::HERMITE:
      // [0, 1] -> [-infty, +infty]
      return M_SQRT2 * (2.0 * xValue - 1.0);
    default:
      return xValue;
  }
}

double OrthogonalBasis1D::evaluate(size_t basisIndex, double xValue) {
#ifdef USE_DAKOTA
  double invNorm = 1. / std::sqrt(basisPoly->norm_squared(basisIndex));
  double normalized_xValue = normalizeInput(xValue);
  return invNorm * basisPoly->type1_value(normalized_xValue, basisIndex);
#else
  std::cerr << "Error in OrthogonalBasis1D::evaluate: "
            << "SG++ was compiled without DAKOTAsupport!\n";
  return false;
#endif
}

} /* namespace combigrid */
} /* namespace sgpp */
