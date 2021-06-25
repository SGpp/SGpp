// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <stdexcept>

#include <sgpp/optimization/fuzzy/FuzzyInterval.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace sgpp {
namespace optimization {

FuzzyInterval::FuzzyInterval(
    double supportLowerBound, double supportUpperBound, size_t numberOfIntegralSamples) :
    supportLowerBound(supportLowerBound),
    supportUpperBound(supportUpperBound),
    numberOfIntegralSamples(numberOfIntegralSamples) {
}

FuzzyInterval::FuzzyInterval(const FuzzyInterval& other) :
    FuzzyInterval(other.supportLowerBound,
                  other.supportUpperBound,
                  other.numberOfIntegralSamples) {
}

FuzzyInterval::~FuzzyInterval() {
}

double FuzzyInterval::computeL1Norm(NormMode normMode) const {
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

      return norm * (ub - lb) / N;
    }
    case NormMode::ViaConfidenceInterval: {
      for (size_t i = 0; i < numberOfIntegralSamples; i++) {
        const double alpha = static_cast<double>(i) / (N - 1.0);
        const double lb = evaluateConfidenceIntervalLowerBound(alpha);
        const double ub = evaluateConfidenceIntervalUpperBound(alpha);
        // (ub - lb) should be >= 0, but just to be sure
        norm += std::abs(ub - lb);
      }

      return norm / N;
    }
  }

  throw std::invalid_argument("Unknown normMode.");
}

double FuzzyInterval::computeL2Norm(NormMode normMode) const {
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

      return std::sqrt((ub - lb) * norm / N);
    }
    case NormMode::ViaConfidenceInterval: {
      for (size_t i = 0; i < numberOfIntegralSamples; i++) {
        const double alpha = static_cast<double>(i) / (N - 1.0);
        const double lb = evaluateConfidenceIntervalLowerBound(alpha);
        const double ub = evaluateConfidenceIntervalUpperBound(alpha);
        norm += (ub - lb) * (ub - lb);
      }

      return std::sqrt(norm / N);
    }
  }

  throw std::invalid_argument("Unknown normMode.");
}

double FuzzyInterval::computeLinfNorm(NormMode normMode) const {
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

      return norm;
    }
    case NormMode::ViaConfidenceInterval: {
      for (size_t i = 0; i < numberOfIntegralSamples; i++) {
        const double alpha = static_cast<double>(i) / (N - 1.0);
        const double lb = evaluateConfidenceIntervalLowerBound(alpha);
        const double ub = evaluateConfidenceIntervalUpperBound(alpha);

        norm = std::max(ub - lb, norm);
      }

      return norm;
    }
  }

  throw std::invalid_argument("Unknown normMode.");
}

double FuzzyInterval::computeL1Error(const FuzzyInterval& other, NormMode normMode) const {
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

      return norm * (ub - lb) / N;
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

      return normLB / N + normUB / N;
    }
  }

  throw std::invalid_argument("Unknown normMode.");
}

double FuzzyInterval::computeL2Error(const FuzzyInterval& other, NormMode normMode) const {
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

      return std::sqrt(norm * (ub - lb) / N);
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

      return std::sqrt(normLB / N) + std::sqrt(normUB / N);
    }
  }

  throw std::invalid_argument("Unknown normMode.");
}

double FuzzyInterval::computeLinfError(const FuzzyInterval& other, NormMode normMode) const {
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
  }

  throw std::invalid_argument("Unknown normMode.");
}

double FuzzyInterval::computeRelativeL1Error(
    const FuzzyInterval& other, NormMode normMode) const {
  const double absoluteError = computeL1Error(other, normMode);
  const double norm = computeL1Norm(normMode);

  return absoluteError / norm;
}

double FuzzyInterval::computeRelativeL2Error(
    const FuzzyInterval& other, NormMode normMode) const {
  const double absoluteError = computeL2Error(other, normMode);
  const double norm = computeL2Norm(normMode);

  return absoluteError / norm;
}

double FuzzyInterval::computeRelativeLinfError(
    const FuzzyInterval& other, NormMode normMode) const {
  const double absoluteError = computeLinfError(other, normMode);
  const double norm = computeLinfNorm(normMode);

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
