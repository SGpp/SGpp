// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/combigrid/functions/AbstractInfiniteFunctionBasis1D.hpp>
#include <sgpp/combigrid/GeneralFunction.hpp>

#include <sgpp/base/tools/json/JSON.hpp>

#ifdef USE_DAKOTA
#include <BasisPolynomial.hpp>
#include <RandomVariable.hpp>
#endif

#include <string>

namespace sgpp {
namespace combigrid {

enum class OrthogonalPolynomialBasisType {
  LEGENDRE,
  JACOBI,
  HERMITE,
  BOUNDED_NORMAL,
  BOUNDED_LOGNORMAL
};

struct OrthogonalPolynomialBasis1DParameters {
  // type
  OrthogonalPolynomialBasisType type_;

  // Beta distribution -> Jacobi polynomials
  double alpha_;
  double beta_;

  // normal distribution -> Hermite polynomials
  double mean_;
  double stddev_;

  // Lognormal distribution -> numerically computed polynomials
  double logmean_;

  // for bounded variables
  double lowerBound_;
  double upperBound_;
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
  explicit OrthogonalPolynomialBasis1D(OrthogonalPolynomialBasis1DConfiguration& config);
  ~OrthogonalPolynomialBasis1D() override;

  double evaluate(size_t basisIndex, double xValue) override;
  double pdf(double xValue);
  double mean();
  double variance();

  double lowerBound();
  double upperBound();

  OrthogonalPolynomialBasis1DConfiguration getConfiguration();

#ifdef USE_DAKOTA
  std::shared_ptr<Pecos::RandomVariable> getRandomVariable();
#endif

  sgpp::combigrid::SingleFunction getWeightFunction();
  size_t numAdditionalQuadraturePoints();

 private:
  double normalizeInput(double xValue);
  OrthogonalPolynomialBasis1DConfiguration config;
#ifdef USE_DAKOTA
  std::shared_ptr<Pecos::BasisPolynomial> basisPoly;
  std::shared_ptr<Pecos::RandomVariable> rv;
#endif
};

} /* namespace combigrid */
} /* namespace sgpp */
