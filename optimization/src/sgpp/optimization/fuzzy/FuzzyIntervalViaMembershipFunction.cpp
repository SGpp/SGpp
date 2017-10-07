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
    double plateauLowerBound, double plateauUpperBound) :
      supportLowerBound(supportLowerBound),
      supportUpperBound(supportUpperBound),
      plateauLowerBound(plateauLowerBound),
      plateauUpperBound(plateauUpperBound) {
}

FuzzyIntervalViaMembershipFunction::~FuzzyIntervalViaMembershipFunction() {
}

double FuzzyIntervalViaMembershipFunction::evaluateConfidenceIntervalLowerBound(
    double alpha) const {
  if (alpha == 0.0) {
    return supportLowerBound;
  } else if (alpha == 1.0) {
    return plateauLowerBound;
  }

  const double tol = 1e-6;
  double xLower = supportLowerBound;
  double xUpper = plateauLowerBound;

  while (xUpper - xLower > tol) {
    const double xCenter = (xLower + xUpper) / 2.0;
    const double alphaCenter = evaluateMembershipFunction(xCenter);

    if (alpha <= alphaCenter) {
      xUpper = xCenter;
    } else {
      xLower = xCenter;
    }
  }

  return (xLower + xUpper) / 2.0;
}

double FuzzyIntervalViaMembershipFunction::evaluateConfidenceIntervalUpperBound(
    double alpha) const {
  if (alpha == 0.0) {
    return supportUpperBound;
  } else if (alpha == 1.0) {
    return plateauUpperBound;
  }

  const double tol = 1e-6;
  double xLower = plateauUpperBound;
  double xUpper = supportUpperBound;

  while (xUpper - xLower > tol) {
    const double xCenter = (xLower + xUpper) / 2.0;
    const double alphaCenter = evaluateMembershipFunction(xCenter);

    if (alpha <= alphaCenter) {
      xLower = xCenter;
    } else {
      xUpper = xCenter;
    }
  }

  return (xLower + xUpper) / 2.0;
}

}  // namespace optimization
}  // namespace sgpp
