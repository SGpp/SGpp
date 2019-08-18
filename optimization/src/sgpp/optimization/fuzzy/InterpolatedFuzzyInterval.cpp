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

double InterpolatedFuzzyInterval::getCoreLowerBound(
    const base::DataVector& xData,
    const base::DataVector& alphaData) {
  const size_t n = xData.getSize();
  double result = std::numeric_limits<double>::quiet_NaN();

  // go through alphaData and determine turning point
  // (first index where alphaData is not monotonically increasing anymore)
  for (size_t i = 0; i < n - 1; i++) {
    if (alphaData[i + 1] < alphaData[i]) {
      // found turning point
      // ==> return last data point where alphaData increased
      return result;
    } else if (alphaData[i + 1] > alphaData[i]) {
      // still monotonically increasing
      result = xData[i + 1];
    }
  }

  // should not happen
  return std::numeric_limits<double>::quiet_NaN();
}

double InterpolatedFuzzyInterval::getCoreUpperBound(
    const base::DataVector& xData,
    const base::DataVector& alphaData) {
  const size_t n = xData.getSize();
  double result = std::numeric_limits<double>::quiet_NaN();

  // go through alphaData from last to first item and determine turning point
  // (last index where alphaData is not monotonically decreasing anymore)
  for (size_t i = n - 1; i > 0; i--) {
    if (alphaData[i - 1] < alphaData[i]) {
      // found turning point
      // ==> return last data point where alphaData decreased
      return result;
    } else if (alphaData[i - 1] > alphaData[i]) {
      // still monotonically decreasing
      result = xData[i - 1];
    }
  }

  // should not happen
  return std::numeric_limits<double>::quiet_NaN();
}

InterpolatedFuzzyInterval* InterpolatedFuzzyInterval::tryDowncast(FuzzyInterval& fuzzyInterval) {
  return dynamic_cast<InterpolatedFuzzyInterval*>(&fuzzyInterval);
}

InterpolatedFuzzyInterval::InterpolatedFuzzyInterval(
    const base::DataVector& xData,
    const base::DataVector& alphaData) :
  FuzzyIntervalViaMembershipFunction(
      xData[0],
      xData[xData.getSize() - 1],
      InterpolatedFuzzyInterval::getCoreLowerBound(xData, alphaData),
      InterpolatedFuzzyInterval::getCoreUpperBound(xData, alphaData)),
  xData(xData),
  alphaData(alphaData) {
}

InterpolatedFuzzyInterval::InterpolatedFuzzyInterval(const InterpolatedFuzzyInterval& other) :
  InterpolatedFuzzyInterval(other.xData, other.alphaData) {
}

InterpolatedFuzzyInterval::~InterpolatedFuzzyInterval() {}

double InterpolatedFuzzyInterval::evaluateMembershipFunction(double x) const {
  // return 0 if x is not in support
  if ((x <= supportLowerBound) || (x >= supportUpperBound)) {
    return 0.0;
  }

  // do a binary search in index space to find the interval [xData[i], xData[i+1]]
  // which contains x (works since x should be sorted)
  size_t i = 0;
  size_t j = xData.getSize() - 1;

  while (j - i > 1) {
    const size_t iCenter = (i + j) / 2;

    if (x < xData[iCenter]) {
      j = iCenter;
    } else {
      i = iCenter;
    }
  }

  if (j == i + 1) {
    // search interval is one interval in xData
    // ==> do linear interpolation between alpha values
    const double t = (x - xData[i]) / (xData[i + 1] - xData[i]);

    return (1.0 - t) * alphaData[i] + t * alphaData[i + 1];
  } else {
    // search interval is only one point
    // ==> return the corresponding alpha value
    return alphaData[i];
  }
}

const base::DataVector& InterpolatedFuzzyInterval::getXData() const {
  return xData;
}

const base::DataVector& InterpolatedFuzzyInterval::getAlphaData() const {
  return alphaData;
}

}  // namespace optimization
}  // namespace sgpp
