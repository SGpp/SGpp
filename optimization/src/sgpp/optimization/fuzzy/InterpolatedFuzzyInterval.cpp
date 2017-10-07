// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/fuzzy/InterpolatedFuzzyInterval.hpp>

#include <algorithm>
#include <limits>

namespace sgpp {
namespace optimization {

double InterpolatedFuzzyInterval::getPlateauLowerBound(
    const base::DataVector& xData,
    const base::DataVector& alphaData) {
  const size_t n = xData.getSize();
  double result = std::numeric_limits<double>::quiet_NaN();

  for (size_t i = 0; i < n - 1; i++) {
    if (alphaData[i + 1] < alphaData[i]) {
      return result;
    } else if (alphaData[i + 1] > alphaData[i]) {
      result = xData[i + 1];
    }
  }

  return std::numeric_limits<double>::quiet_NaN();
}

double InterpolatedFuzzyInterval::getPlateauUpperBound(
    const base::DataVector& xData,
    const base::DataVector& alphaData) {
  const size_t n = xData.getSize();
  double result = std::numeric_limits<double>::quiet_NaN();

  for (size_t i = n - 1; i > 0; i--) {
    if (alphaData[i - 1] < alphaData[i]) {
      return result;
    } else if (alphaData[i - 1] > alphaData[i]) {
      result = xData[i - 1];
    }
  }

  return std::numeric_limits<double>::quiet_NaN();
}

InterpolatedFuzzyInterval::InterpolatedFuzzyInterval(
    const base::DataVector& xData,
    const base::DataVector& alphaData) :
  FuzzyIntervalViaMembershipFunction(
      xData[0],
      xData[xData.getSize() - 1],
      InterpolatedFuzzyInterval::getPlateauLowerBound(xData, alphaData),
      InterpolatedFuzzyInterval::getPlateauUpperBound(xData, alphaData)),
  xData(xData),
  alphaData(alphaData) {
}

InterpolatedFuzzyInterval::~InterpolatedFuzzyInterval() {}

double InterpolatedFuzzyInterval::evaluateMembershipFunction(double x) const {
  const size_t n = xData.getSize();

  if ((x <= xData[0]) || (x >= xData[n - 1])) {
    return 0.0;
  }

  // TODO(valentjn): improve with binary search
  for (size_t i = 0; i < n - 1; i++) {
    if ((xData[i] <= x) && (x <= xData[i + 1])) {
      const double t = (x - xData[i]) / (xData[i + 1] - xData[i]);

      return (1.0 - t) * alphaData[i] + t * alphaData[i + 1];
    }
  }

  return 0.0;
}

const base::DataVector& InterpolatedFuzzyInterval::getXData() const {
  return xData;
}

const base::DataVector& InterpolatedFuzzyInterval::getAlphaData() const {
  return alphaData;
}

}  // namespace optimization
}  // namespace sgpp
