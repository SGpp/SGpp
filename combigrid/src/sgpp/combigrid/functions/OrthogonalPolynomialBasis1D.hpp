// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/combigrid/functions/AbstractInfiniteFunctionBasis1D.hpp>

#ifdef USE_DAKOTA
#include <BasisPolynomial.hpp>
#endif

namespace sgpp {
namespace combigrid {

enum class OrthogonalPolynomialBasisType { LEGENDRE, JACOBI, HERMITE };

class OrthogonalPolynomialBasis1D : public AbstractInfiniteFunctionBasis1D {
 public:
  OrthogonalPolynomialBasis1D();
  explicit OrthogonalPolynomialBasis1D(OrthogonalPolynomialBasisType basisType);
  virtual ~OrthogonalPolynomialBasis1D();

  double evaluate(size_t basisIndex, double xValue) override;

 private:
  double normalizeInput(double xValue);
  OrthogonalPolynomialBasisType basisType;
#ifdef USE_DAKOTA
  std::shared_ptr<Pecos::BasisPolynomial> basisPoly;
#endif
};

} /* namespace combigrid */
} /* namespace sgpp */
