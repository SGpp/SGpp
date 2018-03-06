// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/combigrid/GeneralFunction.hpp>

#ifdef USE_DAKOTA
#include <BasisPolynomial.hpp>
#include <RandomVariable.hpp>
#endif

#include <string>

namespace sgpp {
namespace combigrid {

enum class ProbabilityDensityFunctionType {
  UNIFORM,
  NORMAL,
  BETA,
  BOUNDED_NORMAL,
  BOUNDED_LOGNORMAL
};

struct ProbabilityDensityFunctionParameters {
  // type
  ProbabilityDensityFunctionType type_;

  // Beta distribution
  double alpha_;
  double beta_;

  // normal distribution
  double mean_;
  double stddev_;

  // Lognormal distribution
  double logmean_;

  // for bounded variables
  double lowerBound_;
  double upperBound_;
};

// --------------------------------------------------------------------------
class ProbabilityDensityFunction1D;

class ProbabilityDensityFunction1DConfiguration {
  friend class ProbabilityDensityFunction1D;

 public:
  ProbabilityDensityFunction1DConfiguration();

  //  ProbabilityDensityFunction1DConfiguration* clone();

  void initConfig();

  void setPdfParameters(ProbabilityDensityFunctionParameters newParameters);

  sgpp::combigrid::ProbabilityDensityFunctionParameters pdfParameters;
};

// --------------------------------------------------------------------------

// Wrapper for probability density functions from DAKOTA
class ProbabilityDensityFunction1D {
 public:
  explicit ProbabilityDensityFunction1D(ProbabilityDensityFunction1DConfiguration& config);
  virtual ~ProbabilityDensityFunction1D();

  double normalizeInput(double xValue);
  double pdf(double xValue);
  double mean();
  double variance();

  double lowerBound();
  double upperBound();

  ProbabilityDensityFunction1DConfiguration getConfiguration();

#ifdef USE_DAKOTA
  std::shared_ptr<Pecos::RandomVariable> getRandomVariable();
#endif

  sgpp::combigrid::SingleFunction getWeightFunction();

 private:
  ProbabilityDensityFunction1DConfiguration config;
#ifdef USE_DAKOTA
  std::shared_ptr<Pecos::RandomVariable> rv;
#endif
};

} /* namespace combigrid */
} /* namespace sgpp */
