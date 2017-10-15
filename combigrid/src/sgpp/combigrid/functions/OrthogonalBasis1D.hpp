// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once
#define M_SQRT2 1.41421356237309504880 

#include <sgpp/combigrid/functions/AbstractInfiniteFunctionBasis1D.hpp>

#ifdef USE_DAKOTA
#include <BasisPolynomial.hpp>
#endif

namespace sgpp {
namespace combigrid {

enum class OrthogonalPolynomialBasisType { LEGENDRE, JACOBI, HERMITE };

class OrthogonalBasis1D : public AbstractInfiniteFunctionBasis1D {
 public:
  OrthogonalBasis1D(OrthogonalPolynomialBasisType basisType);
  virtual ~OrthogonalBasis1D();

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
