// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/GeneralFunction.hpp>
#include <sgpp/combigrid/functions/ProbabilityDensityFunction1D.hpp>

#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/tools/json/JSON.hpp>
#include <sgpp/base/tools/json/json_exception.hpp>
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
ProbabilityDensityFunction1DConfiguration::ProbabilityDensityFunction1DConfiguration()
    : json::JSON() {
  initConfig();
}

ProbabilityDensityFunction1DConfiguration::ProbabilityDensityFunction1DConfiguration(
    const std::string& fileName)
    : json::JSON(fileName) {
  initConfig();
  // initialize structs from file
  // configure grid
  try {
    // parameters for jacobi polynomials
    if (this->contains("jacobi_alpha")) pdfParameters.alpha_ = (*this)["jacobi_alpha"].getDouble();
    if (this->contains("jacobi_beta")) pdfParameters.beta_ = (*this)["jacobi_beta"].getDouble();

    // normal distribution
    if (this->contains("normal_mean")) pdfParameters.mean_ = (*this)["normal_mean"].getDouble();
    if (this->contains("normal_stddev"))
      pdfParameters.stddev_ = (*this)["normal_stddev"].getDouble();

    // lognormal distribution
    if (this->contains("lognormal_logmean"))
      pdfParameters.logmean_ = (*this)["lognormal_logmean"].getDouble();

    // bounds
    if (this->contains("lower_bound"))
      pdfParameters.lowerBound_ = (*this)["lower_bound"].getDouble();
    if (this->contains("upper_bound"))
      pdfParameters.upperBound_ = (*this)["upper_bound"].getDouble();
  } catch (json::json_exception& e) {
    std::cout << e.what() << std::endl;
  }
}

void ProbabilityDensityFunction1DConfiguration::initConfig() {
  // set default config
  pdfParameters.type_ = ProbabilityDensityFunctionType::BOUNDED_NORMAL;

  // parameters for Beta distribution
  pdfParameters.alpha_ = 5.0;
  pdfParameters.beta_ = 5.0;

  // parameters for normal distribution
  pdfParameters.mean_ = 0.0;
  pdfParameters.stddev_ = 1.0;

  // parameters for lognormal distribution
  pdfParameters.logmean_ = 0.0;

  // bounds
  pdfParameters.lowerBound_ = 0.0;
  pdfParameters.upperBound_ = 1.0;
}

ProbabilityDensityFunction1DConfiguration* ProbabilityDensityFunction1DConfiguration::clone() {
  ProbabilityDensityFunction1DConfiguration* clone =
      new ProbabilityDensityFunction1DConfiguration(*this);
  return clone;
}

// --------------------------------------------------------------------------------------------

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
    default:
      rv = std::make_shared<Pecos::UniformRandomVariable>(config.pdfParameters.lowerBound_,
                                                          config.pdfParameters.upperBound_);
  }
#endif
}

ProbabilityDensityFunction1D::~ProbabilityDensityFunction1D() {}

double ProbabilityDensityFunction1D::pdf(double xValue) {
#ifdef USE_DAKOTA
  return rv->pdf(xValue);
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
