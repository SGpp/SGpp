// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/fuzzy/FuzzyIntervalViaConfidenceInterval.hpp>

namespace sgpp {
namespace optimization {

FuzzyIntervalViaConfidenceInterval::FuzzyIntervalViaConfidenceInterval(
    double supportLowerBound, double supportUpperBound, size_t numberOfIntegralSamples,
    double binarySearchTolerance) :
        FuzzyInterval(supportLowerBound, supportUpperBound, numberOfIntegralSamples),
        binarySearchTolerance(binarySearchTolerance) {
}

FuzzyIntervalViaConfidenceInterval::FuzzyIntervalViaConfidenceInterval(
    const FuzzyIntervalViaConfidenceInterval& other) :
        FuzzyIntervalViaConfidenceInterval(
            supportLowerBound, supportUpperBound, numberOfIntegralSamples, binarySearchTolerance) {
}

FuzzyIntervalViaConfidenceInterval::~FuzzyIntervalViaConfidenceInterval() {
}

double FuzzyIntervalViaConfidenceInterval::evaluateMembershipFunction(double x) const {
  // determine if x is in the support of the fuzzy interval
  if ((x <= supportLowerBound) || (x >= supportUpperBound)) {
    return 0.0;
  }

  const double coreLowerBound = evaluateConfidenceIntervalLowerBound(1.0);
  const double coreUpperBound = evaluateConfidenceIntervalUpperBound(1.0);

  // determine if x is in the core of the fuzzy interval
  if ((x >= coreLowerBound) && (x <= coreUpperBound)) {
    return 1.0;
  }

  // now we know that x must either be on the left (monotonically increasing) branch
  // or on the right (monotonically decreasing) branch
  const bool isOnLeftBranch = (x < coreLowerBound);

  // do a binary search in alpha space
  double alphaLower = 0.0;
  double alphaUpper = 1.0;

  while (alphaUpper - alphaLower > binarySearchTolerance) {
    const double alphaCenter = (alphaLower + alphaUpper) / 2.0;

    if (isOnLeftBranch) {
      const double xCenter = evaluateConfidenceIntervalLowerBound(alphaCenter);

      if (x <= xCenter) {
        // x is left of the lower bound of the center of the search interval
        // ==> search "lower" alpha half
        alphaUpper = alphaCenter;
      } else {
        // x is right of the lower bound of the center of the search interval
        // ==> search "higher" alpha half
        alphaLower = alphaCenter;
      }
    } else {
      const double xCenter = evaluateConfidenceIntervalUpperBound(alphaCenter);

      if (x <= xCenter) {
        // x is left of the upper bound of the center of the search interval
        // ==> search "higher" alpha half
        alphaLower = alphaCenter;
      } else {
        // x is right of the upper bound of the center of the search interval
        // ==> search "lower" alpha half
        alphaUpper = alphaCenter;
      }
    }
  }

  // return center of last search interval
  return (alphaLower + alphaUpper) / 2.0;
}

double FuzzyIntervalViaConfidenceInterval::getBinarySearchTolerance() const {
  return binarySearchTolerance;
}

void FuzzyIntervalViaConfidenceInterval::setBinarySearchTolerance(double binarySearchTolerance) {
  this->binarySearchTolerance = binarySearchTolerance;
}

}  // namespace optimization
}  // namespace sgpp
