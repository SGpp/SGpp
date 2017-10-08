// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_FUZZY_FUZZYINTERVAL_HPP
#define SGPP_OPTIMIZATION_FUZZY_FUZZYINTERVAL_HPP

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace optimization {

class FuzzyInterval {
 public:
  const size_t DEFAULT_NUMBER_OF_INTEGRAL_SAMPLES = 10000;

  FuzzyInterval();
  virtual ~FuzzyInterval();

  virtual double evaluateMembershipFunction(double x) const = 0;
  virtual double evaluateConfidenceIntervalLowerBound(double alpha) const = 0;
  virtual double evaluateConfidenceIntervalUpperBound(double alpha) const = 0;

  double approximateL1Norm() const;
  double approximateL2Norm() const;
  double approximateLinfNorm() const;
  double approximateL1Error(const FuzzyInterval& other) const;
  double approximateL2Error(const FuzzyInterval& other) const;
  double approximateLinfError(const FuzzyInterval& other) const;
  double approximateRelativeL1Error(const FuzzyInterval& other) const;
  double approximateRelativeL2Error(const FuzzyInterval& other) const;
  double approximateRelativeLinfError(const FuzzyInterval& other) const;

  size_t getNumberOfIntegralSamples() const;
  void setNumberOfIntegralSamples(size_t numberOfIntegralSamples);

 protected:
  size_t numberOfIntegralSamples;
};

}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_FUZZY_FUZZYINTERVAL_HPP */
