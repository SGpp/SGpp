// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_FUZZY_INTERPOLATEDFUZZYINTERVAL_HPP
#define SGPP_OPTIMIZATION_FUZZY_INTERPOLATEDFUZZYINTERVAL_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/fuzzy/FuzzyIntervalViaMembershipFunction.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

namespace sgpp {
namespace optimization {

class InterpolatedFuzzyInterval : public FuzzyIntervalViaMembershipFunction {
 public:
  static double getCoreLowerBound(const base::DataVector& xData,
                                  const base::DataVector& alphaData);
  static double getCoreUpperBound(const base::DataVector& xData,
                                  const base::DataVector& alphaData);

  InterpolatedFuzzyInterval(const base::DataVector& xData, const base::DataVector& alphaData);
  InterpolatedFuzzyInterval(const InterpolatedFuzzyInterval& other);
  ~InterpolatedFuzzyInterval() override;

  double evaluateMembershipFunction(double x) const override;

  const base::DataVector& getXData() const;
  const base::DataVector& getAlphaData() const;

 protected:
  base::DataVector xData;
  base::DataVector alphaData;
};

}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_FUZZY_INTERPOLATEDFUZZYINTERVAL_HPP */

