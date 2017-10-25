// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/functions/OrthogonalPolynomialBasis1D.hpp>
#include <sgpp/base/tools/json/JSON.hpp>
#include <sgpp/base/tools/json/json_exception.hpp>
#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/globaldef.hpp>

#ifdef USE_DAKOTA
#include <HermiteOrthogPolynomial.hpp>
#include <JacobiOrthogPolynomial.hpp>
#include <LegendreOrthogPolynomial.hpp>
#endif

#include <iostream>
#include <string>

#ifdef USE_DAKOTA
#include <BasisPolynomial.hpp>
#include <pecos_global_defs.hpp>
#endif

namespace sgpp {
namespace combigrid {

// --------------------------------------------------------------------------------------------
OrthogonalPolynomialBasis1DConfiguration::OrthogonalPolynomialBasis1DConfiguration()
    : json::JSON() {
  initConfig();
}

OrthogonalPolynomialBasis1DConfiguration::OrthogonalPolynomialBasis1DConfiguration(
    const std::string& fileName)
    : json::JSON(fileName) {
  initConfig();
  // initialize structs from file
  // configure grid
  try {
    // type of the polynomial
    if (this->contains("polynomial_type"))
      polyParameters.type_ = stringToOrthogonalPolynomialType((*this)["polynomial_type"].get());

    // parameters for jacobi polynomials
    if (this->contains("jacobi_alpha")) polyParameters.alpha_ = (*this)["jacobi_alpha"].getDouble();
    if (this->contains("jacobi_beta")) polyParameters.beta_ = (*this)["jacobi_beta"].getDouble();
  } catch (json::json_exception& e) {
    std::cout << e.what() << std::endl;
  }
}

void OrthogonalPolynomialBasis1DConfiguration::initConfig() {
  // set default config
  polyParameters.type_ = OrthogonalPolynomialBasisType::LEGENDRE;

  // parameters for jacobi polynomials
  polyParameters.alpha_ = 5.0;
  polyParameters.beta_ = 5.0;
}

OrthogonalPolynomialBasis1DConfiguration* OrthogonalPolynomialBasis1DConfiguration::clone() {
  OrthogonalPolynomialBasis1DConfiguration* clone =
      new OrthogonalPolynomialBasis1DConfiguration(*this);
  return clone;
}

sgpp::combigrid::OrthogonalPolynomialBasisType
OrthogonalPolynomialBasis1DConfiguration::stringToOrthogonalPolynomialType(
    std::string& polynomialType) {
  if (polynomialType.compare("Legendre") == 0) {
    return sgpp::combigrid::OrthogonalPolynomialBasisType::LEGENDRE;
  } else if (polynomialType.compare("Jacobi") == 0) {
    return sgpp::combigrid::OrthogonalPolynomialBasisType::JACOBI;
  } else if (polynomialType.compare("Hermite") == 0) {
    return sgpp::combigrid::OrthogonalPolynomialBasisType::HERMITE;
  }
  throw sgpp::base::application_exception("polynomial type is unknown");
}

// --------------------------------------------------------------------------------------------

OrthogonalPolynomialBasis1D::OrthogonalPolynomialBasis1D()
    : basisType(OrthogonalPolynomialBasisType::LEGENDRE) {
#ifdef USE_DAKOTA
  basisPoly = std::make_shared<Pecos::BasisPolynomial>(Pecos::LEGENDRE_ORTHOG);
#else
  std::cerr << "Error in OrthogonalBasis1D::evaluate: "
            << "SG++ was compiled without DAKOTAsupport!\n";
#endif
}

OrthogonalPolynomialBasis1D::OrthogonalPolynomialBasis1D(
    OrthogonalPolynomialBasis1DConfiguration& config)
    : basisType(config.polyParameters.type_) {
#ifdef USE_DAKOTA
  switch (basisType) {
    case OrthogonalPolynomialBasisType::HERMITE:
      basisPoly = std::make_shared<Pecos::HermiteOrthogPolynomial>();
      break;
    case OrthogonalPolynomialBasisType::JACOBI:
      basisPoly = std::make_shared<Pecos::JacobiOrthogPolynomial>(config.polyParameters.alpha_,
                                                                  config.polyParameters.beta_);
      break;
    case OrthogonalPolynomialBasisType::LEGENDRE:
      basisPoly = std::make_shared<Pecos::LegendreOrthogPolynomial>();
      break;
    default:
      basisPoly = std::make_shared<Pecos::LegendreOrthogPolynomial>();
  }
#endif
}

OrthogonalPolynomialBasis1D::~OrthogonalPolynomialBasis1D() {}

double OrthogonalPolynomialBasis1D::normalizeInput(double xValue) {
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

double OrthogonalPolynomialBasis1D::evaluate(size_t basisIndex, double xValue) {
#ifdef USE_DAKOTA
  double invNorm = 1. / std::sqrt(basisPoly->norm_squared(static_cast<unsigned short>(basisIndex)));
  double normalized_xValue = normalizeInput(xValue);
  return invNorm *
         basisPoly->type1_value(normalized_xValue, static_cast<unsigned short>(basisIndex));
#else
  std::cerr << "Error in OrthogonalBasis1D::evaluate: "
            << "SG++ was compiled without DAKOTAsupport!\n";
  return false;
#endif
}

} /* namespace combigrid */
} /* namespace sgpp */
