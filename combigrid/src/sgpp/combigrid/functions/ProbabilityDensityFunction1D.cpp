// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/GeneralFunction.hpp>
#include <sgpp/combigrid/functions/ProbabilityDensityFunction1D.hpp>

#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/globaldef.hpp>

#ifdef USE_DAKOTA
#include <pecos_global_defs.hpp>

#include <BetaRandomVariable.hpp>
#include <BoundedLognormalRandomVariable.hpp>
#include <BoundedNormalRandomVariable.hpp>
#include <NormalRandomVariable.hpp>
#include <UniformRandomVariable.hpp>
#endif

#include <iostream>
#include <string>

namespace sgpp {
namespace combigrid {

// --------------------------------------------------------------------------------------------
ProbabilityDensityFunction1DConfiguration::ProbabilityDensityFunction1DConfiguration() {
  initConfig();
}

void ProbabilityDensityFunction1DConfiguration::initConfig() {
  // set default config
  pdfParameters.type_ = ProbabilityDensityFunctionType::BETA;

  // parameters for Beta distribution
  pdfParameters.alpha_ = 3.0;
  pdfParameters.beta_ = 7.0;

  // parameters for normal distribution
  pdfParameters.mean_ = 0.0;
  pdfParameters.stddev_ = 1.0;

  // parameters for lognormal distribution
  pdfParameters.logmean_ = 0.0;

  // bounds
  pdfParameters.lowerBound_ = 0.0;
  pdfParameters.upperBound_ = 1.0;
}

void ProbabilityDensityFunction1DConfiguration::setPdfParameters(
    ProbabilityDensityFunctionParameters newParameters) {
  pdfParameters = newParameters;
}

//  ProbabilityDensityFunction1DConfiguration* ProbabilityDensityFunction1DConfiguration::clone() {
//  ProbabilityDensityFunction1DConfiguration* clone =
//      new ProbabilityDensityFunction1DConfiguration(*this);
//  return clone;
//}

// --------------------------------------------------------------------------------------------

/**
 * @param xValue 	x value in the support of the probability density function
 * @return 			x value transformed to [0,1]
 */
double ProbabilityDensityFunction1D::normalizeInput(double xValue) {
  switch (config.pdfParameters.type_) {
    // UNIFORM, BETA, BOUNDED_NORMAL, NOUNDED_LOGNORMAL, NORMAL require xValue to be in [0,1]
    case ProbabilityDensityFunctionType::UNIFORM:
    case ProbabilityDensityFunctionType::BETA:
    case ProbabilityDensityFunctionType::BOUNDED_NORMAL:
    case ProbabilityDensityFunctionType::BOUNDED_LOGNORMAL:
    case ProbabilityDensityFunctionType::NORMAL:
      break;
  }

  double a = config.pdfParameters.lowerBound_;
  double b = config.pdfParameters.upperBound_;
  if (a == b) {
    std::cerr << "ProbabilityDensityFunction1D: lower and upper bound cannot be identical!"
              << std::endl;
  }
  return a + (b - a) * xValue;
}

ProbabilityDensityFunction1D::ProbabilityDensityFunction1D(
    ProbabilityDensityFunction1DConfiguration& config)
    : config(config) {
#ifdef USE_DAKOTA

  switch (config.pdfParameters.type_) {
    case ProbabilityDensityFunctionType::NORMAL:
      rv = std::make_shared<Pecos::NormalRandomVariable>(config.pdfParameters.mean_,
                                                         config.pdfParameters.stddev_);
      break;
    case ProbabilityDensityFunctionType::BETA:
      rv = std::make_shared<Pecos::BetaRandomVariable>(
          config.pdfParameters.alpha_, config.pdfParameters.beta_, config.pdfParameters.lowerBound_,
          config.pdfParameters.upperBound_);
      break;
    case ProbabilityDensityFunctionType::BOUNDED_LOGNORMAL:
      rv = std::make_shared<Pecos::BoundedLognormalRandomVariable>(
          config.pdfParameters.logmean_, config.pdfParameters.stddev_,
          config.pdfParameters.lowerBound_, config.pdfParameters.upperBound_);
      break;
    case ProbabilityDensityFunctionType::BOUNDED_NORMAL:
      rv = std::make_shared<Pecos::BoundedNormalRandomVariable>(
          config.pdfParameters.mean_, config.pdfParameters.stddev_,
          config.pdfParameters.lowerBound_, config.pdfParameters.upperBound_);
      break;
    case ProbabilityDensityFunctionType::UNIFORM:
    default:
      rv = std::make_shared<Pecos::UniformRandomVariable>(config.pdfParameters.lowerBound_,
                                                          config.pdfParameters.upperBound_);
  }
#endif
}

ProbabilityDensityFunction1D::~ProbabilityDensityFunction1D() {}

double ProbabilityDensityFunction1D::pdf(double xValue) {
#ifdef USE_DAKOTA
  double normalized_xValue = normalizeInput(xValue);
  return rv->pdf(normalized_xValue);
#else
  std::cerr << "Error in ProbabilityDensityFunction1D::pdf: "
            << "SG++ was compiled without DAKOTAsupport!" << std::endl;
  return 0.0;
#endif
}

double ProbabilityDensityFunction1D::mean() {
#ifdef USE_DAKOTA
  return rv->mean();

#else
  std::cerr << "Error in ProbabilityDensityFunction1D::mean: "
            << "SG++ was compiled without DAKOTAsupport!" << std::endl;
  return 0.0;
#endif
}

double ProbabilityDensityFunction1D::variance() {
#ifdef USE_DAKOTA
  return rv->variance();

#else
  std::cerr << "Error in ProbabilityDensityFunction1D::variance: "
            << "SG++ was compiled without DAKOTAsupport!" << std::endl;
  return 0.0;
#endif
}

#ifdef USE_DAKOTA
std::shared_ptr<Pecos::RandomVariable> ProbabilityDensityFunction1D::getRandomVariable() {
  return rv;
}
#endif

sgpp::combigrid::SingleFunction ProbabilityDensityFunction1D::getWeightFunction() {
  return SingleFunction([this](double x_prob) { return this->pdf(x_prob); });
}

double ProbabilityDensityFunction1D::lowerBound() {
#ifdef USE_DAKOTA
  Pecos::RealRealPair bounds_idim = rv->bounds();
  return bounds_idim.first;
#else
  return 0.0;
#endif
}

double ProbabilityDensityFunction1D::upperBound() {
#ifdef USE_DAKOTA
  Pecos::RealRealPair bounds_idim = rv->bounds();
  return bounds_idim.second;
#else
  return 1.0;
#endif
}

ProbabilityDensityFunction1DConfiguration ProbabilityDensityFunction1D::getConfiguration() {
  return config;
}

} /* namespace combigrid */
} /* namespace sgpp */
