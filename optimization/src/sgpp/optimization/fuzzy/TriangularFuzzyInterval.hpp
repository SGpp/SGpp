// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_FUZZY_TRIANGULARFUZZYINTERVAL_HPP
#define SGPP_OPTIMIZATION_FUZZY_TRIANGULARFUZZYINTERVAL_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/fuzzy/FuzzyInterval.hpp>

namespace sgpp {
namespace optimization {

class TriangularFuzzyInterval : public FuzzyInterval {
 public:
  TriangularFuzzyInterval(double mean, double spread);
  TriangularFuzzyInterval(double mean, double leftSpread, double rightSpread);
  TriangularFuzzyInterval(double leftMean, double rightMean,
                          double leftSpread, double rightSpread);

  ~TriangularFuzzyInterval() override;

  double evaluateMembershipFunction(double x) const override;

  double evaluateConfidenceIntervalLowerBound(double alpha) const override;
  double evaluateConfidenceIntervalUpperBound(double alpha) const override;

 protected:
  double leftMean;
  double rightMean;
  double leftSpread;
  double rightSpread;
};

}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_FUZZY_TRIANGULARFUZZYINTERVAL_HPP */
