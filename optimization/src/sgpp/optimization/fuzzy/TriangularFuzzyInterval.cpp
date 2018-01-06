// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/fuzzy/TriangularFuzzyInterval.hpp>

#include <algorithm>

namespace sgpp {
namespace optimization {

TriangularFuzzyInterval::TriangularFuzzyInterval(double mean, double spread) :
  TriangularFuzzyInterval(mean, spread, spread) {
}

TriangularFuzzyInterval::TriangularFuzzyInterval(
    double mean, double leftSpread, double rightSpread) :
  TriangularFuzzyInterval(mean, mean, leftSpread, rightSpread) {
}

TriangularFuzzyInterval::TriangularFuzzyInterval(
    double leftMean, double rightMean, double leftSpread, double rightSpread) :
        FuzzyInterval(leftMean - leftSpread, rightMean + rightSpread),
        leftMean(leftMean), rightMean(rightMean),
        leftSpread(leftSpread), rightSpread(rightSpread) {
}

TriangularFuzzyInterval::TriangularFuzzyInterval(const TriangularFuzzyInterval& other) :
  TriangularFuzzyInterval(other.leftMean, other.rightMean, other.leftSpread, other.rightSpread) {
}

TriangularFuzzyInterval::~TriangularFuzzyInterval() {
}

double TriangularFuzzyInterval::evaluateMembershipFunction(double x) const {
  if (x < leftMean) {
    // x is on the left (monotonically increasing) branch
    return std::max(1.0 - (leftMean - x) / leftSpread, 0.0);
  } else if (x > rightMean) {
    // x is on the right (monotonically dereasing) branch
    return std::max(1.0 - (x - rightMean) / rightSpread, 0.0);
  } else {
    // x is in the core
    return 1.0;
  }
}

double TriangularFuzzyInterval::evaluateConfidenceIntervalLowerBound(double alpha) const {
  return leftMean + (alpha - 1.0) * leftSpread;
}

double TriangularFuzzyInterval::evaluateConfidenceIntervalUpperBound(double alpha) const {
  return rightMean + (1.0 - alpha) * rightSpread;
}

double TriangularFuzzyInterval::getLeftMean() const {
  return leftMean;
}

double TriangularFuzzyInterval::getRightMean() const {
  return rightMean;
}

double TriangularFuzzyInterval::getLeftSpread() const {
  return leftSpread;
}

double TriangularFuzzyInterval::getRightSpread() const {
  return rightSpread;
}

}  // namespace optimization
}  // namespace sgpp
