// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/fuzzy/QuasiGaussianFuzzyNumber.hpp>

#include <algorithm>
#include <cmath>

namespace sgpp {
namespace optimization {

QuasiGaussianFuzzyNumber::QuasiGaussianFuzzyNumber(double mean, double stdev, double cutoff) :
  FuzzyIntervalViaMembershipFunction(mean - cutoff * stdev, mean + cutoff * stdev, mean, mean),
  mean(mean), stdev(stdev), cutoff(cutoff) {
}

QuasiGaussianFuzzyNumber::QuasiGaussianFuzzyNumber(const QuasiGaussianFuzzyNumber& other) :
  QuasiGaussianFuzzyNumber(other.mean, other.stdev, other.cutoff) {
}

QuasiGaussianFuzzyNumber::~QuasiGaussianFuzzyNumber() {
}

double QuasiGaussianFuzzyNumber::evaluateMembershipFunction(double x) const {
  if ((x < supportLowerBound) || (x > supportUpperBound)) {
    // x is not in support
    return 0.0;
  } else {
    // x is in support (not in cutoff area)
    // ==> calculate and return Gauss density
    const double u = std::abs(x - mean) / stdev;

    return std::exp(-u * u / 2.0);
  }
}

double QuasiGaussianFuzzyNumber::getMean() const {
  return mean;
}

double QuasiGaussianFuzzyNumber::getStdev() const {
  return stdev;
}

double QuasiGaussianFuzzyNumber::getCutoff() const {
  return cutoff;
}

}  // namespace optimization
}  // namespace sgpp
