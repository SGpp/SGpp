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

double FuzzyInterval::approximateL1Norm(NormMode normMode) const {
  const double N = static_cast<double>(numberOfIntegralSamples);
  double norm = 0.0;

  switch (normMode) {
    case NormMode::ViaMembershipFunction: {
      const double lb = getSupportLowerBound();
      const double ub = getSupportUpperBound();

      for (size_t i = 0; i < numberOfIntegralSamples; i++) {
        const double x = lb + (ub - lb) * static_cast<double>(i) / (N - 1.0);
        const double mu = evaluateMembershipFunction(x);
        norm += std::abs(mu);
      }

      break;
    }
    case NormMode::ViaConfidenceInterval: {
      for (size_t i = 0; i < numberOfIntegralSamples; i++) {
        const double alpha = static_cast<double>(i) / (N - 1.0);
        const double lb = evaluateConfidenceIntervalLowerBound(alpha);
        const double ub = evaluateConfidenceIntervalUpperBound(alpha);
        // (ub - lb) should be >= 0, but just to be sure
        norm += std::abs(ub - lb);
      }

      break;
    }
    default:
      throw std::invalid_argument("Unknown normMode.");
  }

  norm /= N;

  return norm;
}

double FuzzyInterval::approximateL2Norm(NormMode normMode) const {
  const double N = static_cast<double>(numberOfIntegralSamples);
  double norm = 0.0;

  switch (normMode) {
    case NormMode::ViaMembershipFunction: {
      const double lb = getSupportLowerBound();
      const double ub = getSupportUpperBound();

      for (size_t i = 0; i < numberOfIntegralSamples; i++) {
        const double x = lb + (ub - lb) * static_cast<double>(i) / (N - 1.0);
        const double mu = evaluateMembershipFunction(x);
        norm += mu * mu;
      }

      break;
    }
    case NormMode::ViaConfidenceInterval: {
      for (size_t i = 0; i < numberOfIntegralSamples; i++) {
        const double alpha = static_cast<double>(i) / (N - 1.0);
        const double lb = evaluateConfidenceIntervalLowerBound(alpha);
        const double ub = evaluateConfidenceIntervalUpperBound(alpha);
        norm += (ub - lb) * (ub - lb);
      }

      break;
    }
    default:
      throw std::invalid_argument("Unknown normMode.");
  }

  norm = std::sqrt(norm / N);

  return norm;
}

double FuzzyInterval::approximateLinfNorm(NormMode normMode) const {
  const double N = static_cast<double>(numberOfIntegralSamples);
  double norm = 0.0;

  switch (normMode) {
    case NormMode::ViaMembershipFunction: {
      const double lb = getSupportLowerBound();
      const double ub = getSupportUpperBound();

      for (size_t i = 0; i < numberOfIntegralSamples; i++) {
        const double x = lb + (ub - lb) * static_cast<double>(i) / (N - 1.0);
        const double mu = evaluateMembershipFunction(x);
        norm = std::max(std::abs(mu), norm);
      }

      break;
    }
    case NormMode::ViaConfidenceInterval: {
      for (size_t i = 0; i < numberOfIntegralSamples; i++) {
        const double alpha = static_cast<double>(i) / (N - 1.0);
        const double lb = evaluateConfidenceIntervalLowerBound(alpha);
        const double ub = evaluateConfidenceIntervalUpperBound(alpha);

        norm = std::max(ub - lb, norm);
      }

      break;
    }
    default:
      throw std::invalid_argument("Unknown normMode.");
  }

  return norm;
}

double FuzzyInterval::approximateL1Error(const FuzzyInterval& other, NormMode normMode) const {
  const double N = static_cast<double>(numberOfIntegralSamples);

  switch (normMode) {
    case NormMode::ViaMembershipFunction: {
      const double lb = std::min(getSupportLowerBound(), other.getSupportLowerBound());
      const double ub = std::max(getSupportUpperBound(), other.getSupportUpperBound());
      double norm = 0.0;

      for (size_t i = 0; i < numberOfIntegralSamples; i++) {
        const double x = lb + (ub - lb) * static_cast<double>(i) / (N - 1.0);
        const double mu1 = evaluateMembershipFunction(x);
        const double mu2 = other.evaluateMembershipFunction(x);
        norm += std::abs(mu1 - mu2);
      }

      norm /= N;

      return norm;
    }
    case NormMode::ViaConfidenceInterval: {
      double normLB = 0.0;
      double normUB = 0.0;

      for (size_t i = 0; i < numberOfIntegralSamples; i++) {
        const double alpha = static_cast<double>(i) / (N - 1.0);
        const double lb1 = evaluateConfidenceIntervalLowerBound(alpha);
        const double ub1 = evaluateConfidenceIntervalUpperBound(alpha);
        const double lb2 = other.evaluateConfidenceIntervalLowerBound(alpha);
        const double ub2 = other.evaluateConfidenceIntervalUpperBound(alpha);

        normLB += std::abs(lb1 - lb2);
        normUB += std::abs(ub1 - ub2);
      }

      normLB /= N;
      normUB /= N;

      return normLB + normUB;
    }
    default:
      throw std::invalid_argument("Unknown normMode.");
  }
}

double FuzzyInterval::approximateL2Error(const FuzzyInterval& other, NormMode normMode) const {
  const double N = static_cast<double>(numberOfIntegralSamples);

  switch (normMode) {
    case NormMode::ViaMembershipFunction: {
      const double N = static_cast<double>(numberOfIntegralSamples);
      const double lb = std::min(getSupportLowerBound(), other.getSupportLowerBound());
      const double ub = std::max(getSupportUpperBound(), other.getSupportUpperBound());
      double norm = 0.0;

      for (size_t i = 0; i < numberOfIntegralSamples; i++) {
        const double x = lb + (ub - lb) * static_cast<double>(i) / (N - 1.0);
        const double mu1 = evaluateMembershipFunction(x);
        const double mu2 = other.evaluateMembershipFunction(x);
        norm += (mu1 - mu2) * (mu1 - mu2);
      }

      norm /= N;

      return norm;
    }
    case NormMode::ViaConfidenceInterval: {
      double normLB = 0.0;
      double normUB = 0.0;

      for (size_t i = 0; i < numberOfIntegralSamples; i++) {
        const double alpha = static_cast<double>(i) / (N - 1.0);
        const double lb1 = evaluateConfidenceIntervalLowerBound(alpha);
        const double ub1 = evaluateConfidenceIntervalUpperBound(alpha);
        const double lb2 = other.evaluateConfidenceIntervalLowerBound(alpha);
        const double ub2 = other.evaluateConfidenceIntervalUpperBound(alpha);

        normLB += (lb1 - lb2) * (lb1 - lb2);
        normUB += (ub1 - ub2) * (ub1 - ub2);
      }

      normLB = std::sqrt(normLB / N);
      normUB = std::sqrt(normUB / N);

      return normLB + normUB;
    }
    default:
      throw std::invalid_argument("Unknown normMode.");
  }
}

double FuzzyInterval::approximateLinfError(const FuzzyInterval& other, NormMode normMode) const {
  const double N = static_cast<double>(numberOfIntegralSamples);

  switch (normMode) {
    case NormMode::ViaMembershipFunction: {
      const double lb = std::min(getSupportLowerBound(), other.getSupportLowerBound());
      const double ub = std::max(getSupportUpperBound(), other.getSupportUpperBound());
      double norm = 0.0;

      for (size_t i = 0; i < numberOfIntegralSamples; i++) {
        const double x = lb + (ub - lb) * static_cast<double>(i) / (N - 1.0);
        const double mu1 = evaluateMembershipFunction(x);
        const double mu2 = other.evaluateMembershipFunction(x);
        norm = std::max(std::abs(mu1 - mu2), norm);
      }

      return norm;
    }
    case NormMode::ViaConfidenceInterval: {
      double normLB = 0.0;
      double normUB = 0.0;

      for (size_t i = 0; i < numberOfIntegralSamples; i++) {
        const double alpha = static_cast<double>(i) / (N - 1.0);
        const double lb1 = evaluateConfidenceIntervalLowerBound(alpha);
        const double ub1 = evaluateConfidenceIntervalUpperBound(alpha);
        const double lb2 = other.evaluateConfidenceIntervalLowerBound(alpha);
        const double ub2 = other.evaluateConfidenceIntervalUpperBound(alpha);

        normLB = std::max(std::abs(lb1 - lb2), normLB);
        normUB = std::max(std::abs(ub1 - ub2), normUB);
      }

      return normLB + normUB;
    }
    default:
      throw std::invalid_argument("Unknown normMode.");
  }
}

double FuzzyInterval::approximateRelativeL1Error(
    const FuzzyInterval& other, NormMode normMode) const {
  const double absoluteError = approximateL1Error(other, normMode);
  const double norm = approximateL1Norm(normMode);

  return absoluteError / norm;
}

double FuzzyInterval::approximateRelativeL2Error(
    const FuzzyInterval& other, NormMode normMode) const {
  const double absoluteError = approximateL2Error(other, normMode);
  const double norm = approximateL2Norm(normMode);

  return absoluteError / norm;
}

double FuzzyInterval::approximateRelativeLinfError(
    const FuzzyInterval& other, NormMode normMode) const {
  const double absoluteError = approximateLinfError(other, normMode);
  const double norm = approximateLinfNorm(normMode);

  return absoluteError / norm;
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
