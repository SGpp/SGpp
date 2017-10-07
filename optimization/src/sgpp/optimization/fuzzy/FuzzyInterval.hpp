// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_FUZZY_FUZZYINTERVAL_HPP
#define SGPP_OPTIMIZATION_FUZZY_FUZZYINTERVAL_HPP

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace optimization {

class FuzzyInterval {
 public:
  virtual ~FuzzyInterval() {}

  virtual double evaluateMembershipFunction(double x) const = 0;
  virtual double evaluateConfidenceIntervalLowerBound(double alpha) const = 0;
  virtual double evaluateConfidenceIntervalUpperBound(double alpha) const = 0;
};

}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_FUZZY_FUZZYINTERVAL_HPP */
