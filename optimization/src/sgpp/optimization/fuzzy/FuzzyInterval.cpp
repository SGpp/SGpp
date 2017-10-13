// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/fuzzy/FuzzyInterval.hpp>

#include <algorithm>
#include <cmath>

namespace sgpp {
namespace optimization {

FuzzyInterval::FuzzyInterval(double supportLowerBound, double supportUpperBound) :
    supportLowerBound(supportLowerBound),
    supportUpperBound(supportUpperBound),
    numberOfIntegralSamples(DEFAULT_NUMBER_OF_INTEGRAL_SAMPLES) {
}

FuzzyInterval::~FuzzyInterval() {
}

double FuzzyInterval::approximateL1Norm() const {
  const double N = static_cast<double>(numberOfIntegralSamples);
  double result = 0.0;

  for (size_t i = 0; i < numberOfIntegralSamples; i++) {
    const double alpha = static_cast<double>(i) / (N - 1.0);
    const double lb = evaluateConfidenceIntervalLowerBound(alpha);
    const double ub = evaluateConfidenceIntervalUpperBound(alpha);

    // (ub - lb) should be >= 0, but just to be sure
    result += std::abs(ub - lb);
  }

  return result / N;
}

double FuzzyInterval::approximateL2Norm() const {
  const double N = static_cast<double>(numberOfIntegralSamples);
  double result = 0.0;

  for (size_t i = 0; i < numberOfIntegralSamples; i++) {
    const double alpha = static_cast<double>(i) / (N - 1.0);
    const double lb = evaluateConfidenceIntervalLowerBound(alpha);
    const double ub = evaluateConfidenceIntervalUpperBound(alpha);

    result += (ub - lb) * (ub - lb);
  }

  return std::sqrt(result / N);
}

double FuzzyInterval::approximateLinfNorm() const {
  const double N = static_cast<double>(numberOfIntegralSamples);
  double result = 0.0;

  for (size_t i = 0; i < numberOfIntegralSamples; i++) {
    const double alpha = static_cast<double>(i) / (N - 1.0);
    const double lb = evaluateConfidenceIntervalLowerBound(alpha);
    const double ub = evaluateConfidenceIntervalUpperBound(alpha);

    result = std::max(ub - lb, result);
  }

  return result;
}

double FuzzyInterval::approximateL1Error(const FuzzyInterval& other) const {
  const double N = static_cast<double>(numberOfIntegralSamples);
  double result = 0.0;

  for (size_t i = 0; i < numberOfIntegralSamples; i++) {
    const double alpha = static_cast<double>(i) / (N - 1.0);
    const double lb1 = evaluateConfidenceIntervalLowerBound(alpha);
    const double ub1 = evaluateConfidenceIntervalUpperBound(alpha);
    const double lb2 = other.evaluateConfidenceIntervalLowerBound(alpha);
    const double ub2 = other.evaluateConfidenceIntervalUpperBound(alpha);

    result += std::abs(lb1 - lb2) + std::abs(ub1 - ub2);
  }

  return result / N;
}

double FuzzyInterval::approximateL2Error(const FuzzyInterval& other) const {
  const double N = static_cast<double>(numberOfIntegralSamples);
  double result = 0.0;

  for (size_t i = 0; i < numberOfIntegralSamples; i++) {
    const double alpha = static_cast<double>(i) / (N - 1.0);
    const double lb1 = evaluateConfidenceIntervalLowerBound(alpha);
    const double ub1 = evaluateConfidenceIntervalUpperBound(alpha);
    const double lb2 = other.evaluateConfidenceIntervalLowerBound(alpha);
    const double ub2 = other.evaluateConfidenceIntervalUpperBound(alpha);

    result += (lb1 - lb2) * (lb1 - lb2) + (ub1 - ub2) * (ub1 - ub2);
  }

  return std::sqrt(result / N);
}

double FuzzyInterval::approximateLinfError(const FuzzyInterval& other) const {
  const double N = static_cast<double>(numberOfIntegralSamples);
  double resultLB = 0.0;
  double resultUB = 0.0;

  for (size_t i = 0; i < numberOfIntegralSamples; i++) {
    const double alpha = static_cast<double>(i) / (N - 1.0);
    const double lb1 = evaluateConfidenceIntervalLowerBound(alpha);
    const double ub1 = evaluateConfidenceIntervalUpperBound(alpha);
    const double lb2 = other.evaluateConfidenceIntervalLowerBound(alpha);
    const double ub2 = other.evaluateConfidenceIntervalUpperBound(alpha);

    resultLB = std::max(std::abs(lb1 - lb2), resultLB);
    resultUB = std::max(std::abs(ub1 - ub2), resultUB);
  }

  return std::max(resultLB, resultUB);
}

double FuzzyInterval::approximateRelativeL1Error(const FuzzyInterval& other) const {
  const double N = static_cast<double>(numberOfIntegralSamples);
  double normDifference = 0.0;
  double normAbsolute = 0.0;

  for (size_t i = 0; i < numberOfIntegralSamples; i++) {
    const double alpha = static_cast<double>(i) / (N - 1.0);
    const double lb1 = evaluateConfidenceIntervalLowerBound(alpha);
    const double ub1 = evaluateConfidenceIntervalUpperBound(alpha);
    const double lb2 = other.evaluateConfidenceIntervalLowerBound(alpha);
    const double ub2 = other.evaluateConfidenceIntervalUpperBound(alpha);

    normDifference += std::abs(lb1 - lb2) + std::abs(ub1 - ub2);
    // (ub1 - lb1) should be >= 0, but just to be sure
    normAbsolute += std::abs(ub1 - lb1);
  }

  return normDifference / normAbsolute;
}

double FuzzyInterval::approximateRelativeL2Error(const FuzzyInterval& other) const {
  const double N = static_cast<double>(numberOfIntegralSamples);
  double normDifference = 0.0;
  double normAbsolute = 0.0;

  for (size_t i = 0; i < numberOfIntegralSamples; i++) {
    const double alpha = static_cast<double>(i) / (N - 1.0);
    const double lb1 = evaluateConfidenceIntervalLowerBound(alpha);
    const double ub1 = evaluateConfidenceIntervalUpperBound(alpha);
    const double lb2 = other.evaluateConfidenceIntervalLowerBound(alpha);
    const double ub2 = other.evaluateConfidenceIntervalUpperBound(alpha);

    normDifference += (lb1 - lb2) * (lb1 - lb2) + (ub1 - ub2) * (ub1 - ub2);
    normAbsolute += (ub1 - lb1) * (ub1 - lb1);
  }

  return std::sqrt(normDifference / normAbsolute);
}

double FuzzyInterval::approximateRelativeLinfError(const FuzzyInterval& other) const {
  const double N = static_cast<double>(numberOfIntegralSamples);
  double normDifferenceLB = 0.0;
  double normDifferenceUB = 0.0;
  double normAbsolute = 0.0;

  for (size_t i = 0; i < numberOfIntegralSamples; i++) {
    const double alpha = static_cast<double>(i) / (N - 1.0);
    const double lb1 = evaluateConfidenceIntervalLowerBound(alpha);
    const double ub1 = evaluateConfidenceIntervalUpperBound(alpha);
    const double lb2 = other.evaluateConfidenceIntervalLowerBound(alpha);
    const double ub2 = other.evaluateConfidenceIntervalUpperBound(alpha);

    normDifferenceLB = std::max(std::abs(lb1 - lb2), normDifferenceLB);
    normDifferenceUB = std::max(std::abs(ub1 - ub2), normDifferenceUB);
    normAbsolute = std::max(ub1 - lb1, normAbsolute);
  }

  return std::max(normDifferenceLB, normDifferenceUB) / normAbsolute;
}

double FuzzyInterval::getSupportLowerBound() const {
  return supportLowerBound;
}

double FuzzyInterval::getSupportUpperBound() const {
  return supportUpperBound;
}

size_t FuzzyInterval::getNumberOfIntegralSamples() const {
  return numberOfIntegralSamples;
}

void FuzzyInterval::setNumberOfIntegralSamples(size_t numberOfIntegralSamples) {
  this->numberOfIntegralSamples = numberOfIntegralSamples;
}

}  // namespace optimization
}  // namespace sgpp
