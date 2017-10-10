// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_FUZZY_FUZZYINTERVALVIACONFIDENCEINTERVAL_HPP
#define SGPP_OPTIMIZATION_FUZZY_FUZZYINTERVALVIACONFIDENCEINTERVAL_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/fuzzy/FuzzyInterval.hpp>

namespace sgpp {
namespace optimization {

class FuzzyIntervalViaConfidenceInterval : public FuzzyInterval {
 public:
  FuzzyIntervalViaConfidenceInterval(double supportLowerBound, double supportUpperBound);
  ~FuzzyIntervalViaConfidenceInterval() override;

  double evaluateMembershipFunction(double x) const override;
};

}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_FUZZY_FUZZYINTERVALVIACONFIDENCEINTERVAL_HPP */

