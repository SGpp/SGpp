// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/functions/OrthogonalPolynomialBasis1D.hpp>
#include <sgpp/combigrid/GeneralFunction.hpp>

#include <sgpp/base/tools/json/JSON.hpp>
#include <sgpp/base/tools/json/json_exception.hpp>
#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/globaldef.hpp>

#ifdef USE_DAKOTA
#include <pecos_global_defs.hpp>
#include <BasisPolynomial.hpp>
#include <HermiteOrthogPolynomial.hpp>
#include <JacobiOrthogPolynomial.hpp>
#include <LegendreOrthogPolynomial.hpp>
#include <NumericGenOrthogPolynomial.hpp>

#include <NormalRandomVariable.hpp>
#include <UniformRandomVariable.hpp>
#include <BetaRandomVariable.hpp>
#include <BoundedLognormalRandomVariable.hpp>
#endif

#include <string>

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

    // parameters for orthogonal polynomials wrt normal distribution
    if (this->contains("normal_mean")) polyParameters.mean_ = (*this)["normal_mean"].getDouble();
    if (this->contains("normal_stddev"))
      polyParameters.stddev_ = (*this)["normal_stddev"].getDouble();

    // parameters for orthogonal polynomials wrt lognormal distribution
    if (this->contains("lognormal_logmean"))
      polyParameters.logmean_ = (*this)["lognormal_logmean"].getDouble();

    // bounds
    if (this->contains("lower_bound"))
      polyParameters.lowerBound_ = (*this)["lower_bound"].getDouble();
    if (this->contains("upper_bound"))
      polyParameters.upperBound_ = (*this)["upper_bound"].getDouble();
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

  // parameters for normal distribution
  polyParameters.mean_ = 0.0;
  polyParameters.stddev_ = 1.0;

  // parameters for lognormal distribution
  polyParameters.logmean_ = 0.0;

  // bounds
  polyParameters.lowerBound_ = 0.0;
  polyParameters.upperBound_ = 1.0;
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
  } else if (polynomialType.compare("Bounded_Lognormal") == 0) {
    return sgpp::combigrid::OrthogonalPolynomialBasisType::BOUNDED_LOGNORMAL;
  } else if (polynomialType.compare("Bounded_Normal") == 0) {
    return sgpp::combigrid::OrthogonalPolynomialBasisType::BOUNDED_NORMAL;
  } else {
    throw sgpp::base::application_exception("polynomial type is unknown");
  }
}

// --------------------------------------------------------------------------------------------

OrthogonalPolynomialBasis1D::OrthogonalPolynomialBasis1D(
    OrthogonalPolynomialBasis1DConfiguration& config)
    : config(config) {
#ifdef USE_DAKOTA
  // pointers are not owned by the objects, so do not delete them
  Pecos::NumericGenOrthogPolynomial* numericPoly;

  switch (config.polyParameters.type_) {
    case OrthogonalPolynomialBasisType::HERMITE:
      basisPoly = std::make_shared<Pecos::BasisPolynomial>(Pecos::HERMITE_ORTHOG);
      rv = std::make_shared<Pecos::NormalRandomVariable>(config.polyParameters.mean_,
                                                         config.polyParameters.stddev_);
      break;
    case OrthogonalPolynomialBasisType::JACOBI:
      basisPoly = std::make_shared<Pecos::BasisPolynomial>(Pecos::JACOBI_ORTHOG);
      // set parameters
      basisPoly->alpha_stat(config.polyParameters.alpha_);
      basisPoly->beta_stat(config.polyParameters.beta_);

      rv = std::make_shared<Pecos::BetaRandomVariable>(
          config.polyParameters.alpha_, config.polyParameters.beta_,
          config.polyParameters.lowerBound_, config.polyParameters.upperBound_);
      break;
    case OrthogonalPolynomialBasisType::BOUNDED_LOGNORMAL:
      basisPoly = std::make_shared<Pecos::BasisPolynomial>(Pecos::NUM_GEN_ORTHOG);
      numericPoly = dynamic_cast<Pecos::NumericGenOrthogPolynomial*>(basisPoly->polynomial_rep());
      numericPoly->bounded_lognormal_distribution(
          config.polyParameters.logmean_, config.polyParameters.stddev_,
          config.polyParameters.lowerBound_, config.polyParameters.upperBound_);
      numericPoly->coefficients_norms_flag(true);

      rv = std::make_shared<Pecos::BoundedLognormalRandomVariable>(
          config.polyParameters.logmean_, config.polyParameters.stddev_,
          config.polyParameters.lowerBound_, config.polyParameters.upperBound_);
      break;
    case OrthogonalPolynomialBasisType::BOUNDED_NORMAL:
      basisPoly = std::make_shared<Pecos::BasisPolynomial>(Pecos::HERMITE_ORTHOG);
      rv = std::make_shared<Pecos::BoundedNormalRandomVariable>(
          config.polyParameters.mean_, config.polyParameters.stddev_,
          config.polyParameters.lowerBound_, config.polyParameters.upperBound_);
      break;
    case OrthogonalPolynomialBasisType::LEGENDRE:
    default:
      basisPoly = std::make_shared<Pecos::BasisPolynomial>(Pecos::LEGENDRE_ORTHOG);
      rv = std::make_shared<Pecos::UniformRandomVariable>(config.polyParameters.lowerBound_,
                                                          config.polyParameters.upperBound_);
  }
#endif
}

OrthogonalPolynomialBasis1D::~OrthogonalPolynomialBasis1D() {}

double OrthogonalPolynomialBasis1D::normalizeInput(double xValue) {
  switch (config.polyParameters.type_) {
    case OrthogonalPolynomialBasisType::LEGENDRE:
    case OrthogonalPolynomialBasisType::JACOBI:
      // [0, 1] -> [-1, 1]
      return 2.0 * xValue - 1.0;
    case OrthogonalPolynomialBasisType::BOUNDED_NORMAL:
    case OrthogonalPolynomialBasisType::BOUNDED_LOGNORMAL:
      return xValue;
    //      // [0, 1] -> [lower bound, upper bound]
    //      return (config.polyParameters.upperBound_ - config.polyParameters.lowerBound_) * xValue
    //      +
    //             config.polyParameters.lowerBound_;
    case OrthogonalPolynomialBasisType::HERMITE:
      // [0, 1] -> [-1, 1] -> [-infty, +infty]
      return M_SQRT2 * (2.0 * xValue - 1.0);
  }
}

double OrthogonalPolynomialBasis1D::evaluate(size_t basisIndex, double xValue) {
#ifdef USE_DAKOTA
  uint16_t castedBasisIndex = static_cast<uint16_t>(basisIndex);
  double invNorm = 1. / std::sqrt(basisPoly->norm_squared(castedBasisIndex));
  double normalized_xValue = normalizeInput(xValue);
  return invNorm * basisPoly->type1_value(normalized_xValue, castedBasisIndex);
#else
  throw sgpp::base::application_exception(
      "Error in OrthogonalBasis1D::evaluate: SG++ was compiled without DAKOTA support!");
#endif
}

double OrthogonalPolynomialBasis1D::pdf(double xValue) {
#ifdef USE_DAKOTA
  return rv->pdf(xValue);
#else
  throw sgpp::base::application_exception(
      "Error in OrthogonalBasis1D::evaluate: SG++ was compiled without DAKOTA support!");
#endif
}

double OrthogonalPolynomialBasis1D::mean() {
#ifdef USE_DAKOTA
  return rv->mean();

#else
  throw sgpp::base::application_exception(
      "Error in OrthogonalBasis1D::evaluate: SG++ was compiled without DAKOTA support!");
#endif
}

double OrthogonalPolynomialBasis1D::variance() {
#ifdef USE_DAKOTA
  return rv->variance();

#else
  throw sgpp::base::application_exception(
      "Error in OrthogonalBasis1D::evaluate: SG++ was compiled without DAKOTA support!");
#endif
}

#ifdef USE_DAKOTA
std::shared_ptr<Pecos::RandomVariable> OrthogonalPolynomialBasis1D::getRandomVariable() {
  return rv;
}
#endif

sgpp::combigrid::SingleFunction OrthogonalPolynomialBasis1D::getWeightFunction() {
  return SingleFunction([this](double x_prob) { return this->pdf(x_prob); });
}

size_t OrthogonalPolynomialBasis1D::numAdditionalQuadraturePoints() {
  switch (config.polyParameters.type_) {
    case OrthogonalPolynomialBasisType::LEGENDRE:
      return 5;
    case OrthogonalPolynomialBasisType::BOUNDED_LOGNORMAL:
      return 15;
    case OrthogonalPolynomialBasisType::JACOBI:
      return 5;
    case OrthogonalPolynomialBasisType::HERMITE:
    case OrthogonalPolynomialBasisType::BOUNDED_NORMAL:
      return 10;
  }
}

double OrthogonalPolynomialBasis1D::lowerBound() {
#ifdef USE_DAKOTA
  Pecos::RealRealPair bounds_idim = rv->bounds();
  return bounds_idim.first;
#else
  return 0.0;
#endif
}

double OrthogonalPolynomialBasis1D::upperBound() {
#ifdef USE_DAKOTA
  Pecos::RealRealPair bounds_idim = rv->bounds();
  return bounds_idim.second;
#else
  return 1.0;
#endif
}

OrthogonalPolynomialBasis1DConfiguration OrthogonalPolynomialBasis1D::getConfiguration() {
  return config;
}

} /* namespace combigrid */
} /* namespace sgpp */
