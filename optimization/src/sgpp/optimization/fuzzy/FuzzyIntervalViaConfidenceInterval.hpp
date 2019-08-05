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

/**
 * Abstract class for a fuzzy interval which is defined by stating its
 * confidence intervals \f$(\tilde{x})_\alpha\f$ for all \f$\alpha \in [0, 1]\f$.
 */
class FuzzyIntervalViaConfidenceInterval : public FuzzyInterval {
 public:
  /// default tolerance for the binary search
  static constexpr double DEFAULT_BINARY_SEARCH_TOLERANCE = 1e-6;

  /**
   * Constructor. Needs the support of the fuzzy interval (which is
   * always a closed interval, so it suffices to supply lower and upper bound).
   *
   * @param supportLowerBound       lower bound of the support
   * @param supportUpperBound       upper bound of the support
   * @param binarySearchTolerance   tolerance for the binary search
   */
  FuzzyIntervalViaConfidenceInterval(
      double supportLowerBound, double supportUpperBound,
      double binarySearchTolerance = DEFAULT_BINARY_SEARCH_TOLERANCE);

  /**
   * Destructor.
   */
  ~FuzzyIntervalViaConfidenceInterval() override;

  /**
   * Evaluate the membership function.
   *
   * @param x   \f$x \in X\f$
   * @return    \f$\mu_{\tilde{x}}(x) \in [0, 1]\f$
   */
  double evaluateMembershipFunction(double x) const override;

  /**
   * @return  tolerance for the binary search
   */
  double getBinarySearchTolerance() const;

  /**
   * @param binarySearchTolerance   tolerance for the binary search
   */
  void setBinarySearchTolerance(double binarySearchTolerance);

 protected:
  /// tolerance for the binary search
  double binarySearchTolerance;
};

}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_FUZZY_FUZZYINTERVALVIACONFIDENCEINTERVAL_HPP */

