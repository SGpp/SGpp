// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/fuzzy/FuzzyIntervalViaConfidenceInterval.hpp>

namespace sgpp {
namespace optimization {

FuzzyIntervalViaConfidenceInterval::~FuzzyIntervalViaConfidenceInterval() {
}

double FuzzyIntervalViaConfidenceInterval::evaluateMembershipFunction(double x) const {
  const double supportLowerBound = evaluateConfidenceIntervalLowerBound(0.0);
  const double supportUpperBound = evaluateConfidenceIntervalUpperBound(0.0);

  if ((x <= supportLowerBound) || (x >= supportUpperBound)) {
    return 0.0;
  }

  const double plateauLowerBound = evaluateConfidenceIntervalLowerBound(1.0);
  const double plateauUpperBound = evaluateConfidenceIntervalUpperBound(1.0);

  if ((x >= plateauLowerBound) && (x <= plateauUpperBound)) {
    return 1.0;
  }

  const bool isOnLeftBranch = (x < plateauLowerBound);
  const double tol = 1e-6;
  double alphaLower = 0.0;
  double alphaUpper = 1.0;

  while (alphaUpper - alphaLower > tol) {
    const double alphaCenter = (alphaLower + alphaUpper) / 2.0;

    if (isOnLeftBranch) {
      const double xCenter = evaluateConfidenceIntervalLowerBound(alphaCenter);

      if (x <= xCenter) {
        alphaUpper = alphaCenter;
      } else {
        alphaLower = alphaCenter;
      }
    } else {
      const double xCenter = evaluateConfidenceIntervalUpperBound(alphaCenter);

      if (x <= xCenter) {
        alphaLower = alphaCenter;
      } else {
        alphaUpper = alphaCenter;
      }
    }
  }

  return (alphaLower + alphaUpper) / 2.0;
}

}  // namespace optimization
}  // namespace sgpp
