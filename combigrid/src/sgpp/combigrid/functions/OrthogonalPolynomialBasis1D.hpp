// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/combigrid/functions/AbstractInfiniteFunctionBasis1D.hpp>
#include <sgpp/base/tools/json/JSON.hpp>

#ifdef USE_DAKOTA
#include <BasisPolynomial.hpp>
#endif

namespace sgpp {
namespace combigrid {

enum class OrthogonalPolynomialBasisType { LEGENDRE, JACOBI, HERMITE };

struct OrthogonalPolynomialBasis1DParameters {
  // type
  OrthogonalPolynomialBasisType type_;

  // Jacobi polynomials
  double alpha_;
  double beta_;
};

// --------------------------------------------------------------------------
class OrthogonalPolynomialBasis1D;

class OrthogonalPolynomialBasis1DConfiguration : public json::JSON {
  friend class OrthogonalPolynomialBasis1D;

 public:
  OrthogonalPolynomialBasis1DConfiguration();
  explicit OrthogonalPolynomialBasis1DConfiguration(const std::string& fileName);

  OrthogonalPolynomialBasis1DConfiguration* clone() override;

  void initConfig();
  sgpp::combigrid::OrthogonalPolynomialBasisType stringToOrthogonalPolynomialType(
      std::string& polynomialType);

  sgpp::combigrid::OrthogonalPolynomialBasis1DParameters polyParameters;
};

// --------------------------------------------------------------------------
class OrthogonalPolynomialBasis1D : public AbstractInfiniteFunctionBasis1D {
 public:
  OrthogonalPolynomialBasis1D();
  explicit OrthogonalPolynomialBasis1D(OrthogonalPolynomialBasis1DConfiguration& config);
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
