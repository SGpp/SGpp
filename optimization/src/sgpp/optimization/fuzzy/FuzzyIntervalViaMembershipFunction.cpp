// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/fuzzy/FuzzyIntervalViaMembershipFunction.hpp>

namespace sgpp {
namespace optimization {

FuzzyIntervalViaMembershipFunction::FuzzyIntervalViaMembershipFunction(
    double supportLowerBound, double supportUpperBound,
    double coreLowerBound, double coreUpperBound) :
      supportLowerBound(supportLowerBound),
      supportUpperBound(supportUpperBound),
      coreLowerBound(coreLowerBound),
      coreUpperBound(coreUpperBound) {
}

FuzzyIntervalViaMembershipFunction::~FuzzyIntervalViaMembershipFunction() {
}

double FuzzyIntervalViaMembershipFunction::evaluateConfidenceIntervalLowerBound(
    double alpha) const {
  // confidence interval for alpha = 0 or 1 is known
  if (alpha == 0.0) {
    return supportLowerBound;
  } else if (alpha == 1.0) {
    return coreLowerBound;
  }

  // do a binary search in x space
  const double tol = 1e-6;
  double xLower = supportLowerBound;
  double xUpper = coreLowerBound;

  while (xUpper - xLower > tol) {
    const double xCenter = (xLower + xUpper) / 2.0;
    const double alphaCenter = evaluateMembershipFunction(xCenter);

    if (alpha <= alphaCenter) {
      // alpha is left of alpha value corresponding to the center of the search interval
      // ==> search "lower" x half
      xUpper = xCenter;
    } else {
      // alpha is right of alpha value corresponding to the center of the search interval
      // ==> search "higher" x half
      xLower = xCenter;
    }
  }

  // return center of last search interval
  return (xLower + xUpper) / 2.0;
}

double FuzzyIntervalViaMembershipFunction::evaluateConfidenceIntervalUpperBound(
    double alpha) const {
  // confidence interval for alpha = 0 or 1 is known
  if (alpha == 0.0) {
    return supportUpperBound;
  } else if (alpha == 1.0) {
    return coreUpperBound;
  }

  // do a binary search in x space
  const double tol = 1e-6;
  double xLower = coreUpperBound;
  double xUpper = supportUpperBound;

  while (xUpper - xLower > tol) {
    const double xCenter = (xLower + xUpper) / 2.0;
    const double alphaCenter = evaluateMembershipFunction(xCenter);

    if (alpha <= alphaCenter) {
      // alpha is left of alpha value corresponding to the center of the search interval
      // ==> search "higher" x half
      xLower = xCenter;
    } else {
      // alpha is right of alpha value corresponding to the center of the search interval
      // ==> search "lower" x half
      xUpper = xCenter;
    }
  }

  // return center of last search interval
  return (xLower + xUpper) / 2.0;
}

}  // namespace optimization
}  // namespace sgpp
