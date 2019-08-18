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

/**
 * Abstract class for a fuzzy interval which is defined by stating its
 * membership function \f$\mu_{\tilde{x}}\colon X \to [0, 1]\f$.
 */
class FuzzyIntervalViaMembershipFunction : public FuzzyInterval {
 public:
  /// default tolerance for the binary search
  static constexpr double DEFAULT_BINARY_SEARCH_TOLERANCE = 1e-6;

  /**
   * Constructor. Needs the support of the fuzzy interval and its core, i.e.,
   * \f$(\tilde{x})_\alpha\f$ for \f$alpha = 0\f$ and for \f$\alpha = 1\f$,
   * (which are always closed intervals, so it suffices to supply lower and upper bounds).
   * The core is needed since the binary search cannot be performed for \f$\alpha = 1\f$.
   *
   * @param supportLowerBound         lower bound of the support
   * @param supportUpperBound         upper bound of the support
   * @param coreLowerBound            lower bound of the core
   * @param coreUpperBound            upper bound of the core
   * @param numberOfIntegralSamples   number of samples to compute norms
   * @param binarySearchTolerance     tolerance for the binary search
   */
  FuzzyIntervalViaMembershipFunction(
      double supportLowerBound, double supportUpperBound,
      double coreLowerBound, double coreUpperBound,
      size_t numberOfIntegralSamples = DEFAULT_NUMBER_OF_INTEGRAL_SAMPLES,
      double binarySearchTolerance = DEFAULT_BINARY_SEARCH_TOLERANCE);

  /**
   * Copy constructor.
   *
   * @param other   other fuzzy interval
   */
  FuzzyIntervalViaMembershipFunction(const FuzzyIntervalViaMembershipFunction& other);

  /**
   * Destructor.
   */
  ~FuzzyIntervalViaMembershipFunction() override;

  /**
   * Evaluate the lower bound of a confidence interval,
   * which is always a closed interval \f$(\tilde{x})_\alpha = [a, b]\f$.
   *
   * @param alpha   \f$\alpha \in [0, 1]\f$
   * @return        \f$a \in X\f$
   */
  double evaluateConfidenceIntervalLowerBound(double alpha) const override;

  /**
   * Evaluate the upper bound of a confidence interval,
   * which is always a closed interval \f$(\tilde{x})_\alpha = [a, b]\f$.
   *
   * @param alpha   \f$\alpha \in [0, 1]\f$
   * @return        \f$b \in X\f$
   */
  double evaluateConfidenceIntervalUpperBound(double alpha) const override;

  /**
   * @return  lower bound of the core
   */
  double getCoreLowerBound() const;

  /**
   * @return  upper bound of the core
   */
  double getCoreUpperBound() const;

  /**
   * @return  tolerance for the binary search
   */
  double getBinarySearchTolerance() const;

  /**
   * @param binarySearchTolerance   tolerance for the binary search
   */
  void setBinarySearchTolerance(double binarySearchTolerance);

 protected:
  /// lower bound of the core
  double coreLowerBound;
  /// upper bound of the core
  double coreUpperBound;
  /// tolerance for the binary search
  double binarySearchTolerance;
};

}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_FUZZY_FUZZYINTERVALVIAMEMBERSHIPFUNCTION_HPP */
