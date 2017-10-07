// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_FUZZY_QUASIGAUSSIANFUZZYNUMBER_HPP
#define SGPP_OPTIMIZATION_FUZZY_QUASIGAUSSIANFUZZYNUMBER_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/fuzzy/FuzzyIntervalViaMembershipFunction.hpp>

namespace sgpp {
namespace optimization {

class QuasiGaussianFuzzyNumber : public FuzzyIntervalViaMembershipFunction {
 public:
  QuasiGaussianFuzzyNumber(double mean, double stdev, double cutoff);
  ~QuasiGaussianFuzzyNumber() override;

  double evaluateMembershipFunction(double x) const override;

 protected:
  double mean;
  double stdev;
  double cutoff;
};

}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_FUZZY_QUASIGAUSSIANFUZZYNUMBER_HPP */
