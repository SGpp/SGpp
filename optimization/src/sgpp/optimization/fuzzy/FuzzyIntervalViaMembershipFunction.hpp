// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_FUZZY_FUZZYINTERVALVIAMEMBERSHIPFUNCTION_HPP
#define SGPP_OPTIMIZATION_FUZZY_FUZZYINTERVALVIAMEMBERSHIPFUNCTION_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/fuzzy/FuzzyInterval.hpp>

namespace sgpp {
namespace optimization {

class FuzzyIntervalViaMembershipFunction : public FuzzyInterval {
 public:
  FuzzyIntervalViaMembershipFunction(
      double supportLowerBound, double supportUpperBound,
      double coreLowerBound, double coreUpperBound);
  ~FuzzyIntervalViaMembershipFunction() override;

  double evaluateConfidenceIntervalLowerBound(double alpha) const override;
  double evaluateConfidenceIntervalUpperBound(double alpha) const override;

 protected:
  double supportLowerBound;
  double supportUpperBound;
  double coreLowerBound;
  double coreUpperBound;
};

}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_FUZZY_FUZZYINTERVALVIAMEMBERSHIPFUNCTION_HPP */

