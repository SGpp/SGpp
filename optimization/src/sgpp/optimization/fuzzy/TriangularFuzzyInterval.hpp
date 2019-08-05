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

/**
 * Triangular fuzzy interval; its membership function linearly increases from 0
 * to 1, stays 1, and linearly decreases back to 0.
 * The core (i.e., the area where the membership function equals 1) is given by
 * \f$[\mathrm{leftMean}, \mathrm{rightMean}]\f$.
 * The support is given by
 * \f$[\mathrm{leftMean} - \mathrm{leftSpread}, \mathrm{rightMean} + \mathrm{rightSpread}]\f$.
 */
class TriangularFuzzyInterval : public FuzzyInterval {
 public:
  /**
   * Constructor.
   *
   * @param mean    left mean = right mean
   * @param spread  left spread = right spread
   */
  TriangularFuzzyInterval(double mean, double spread);

  /**
   * Constructor.
   *
   * @param mean          left mean = right mean
   * @param leftSpread    left spread
   * @param rightSpread   right spread
   */
  TriangularFuzzyInterval(double mean, double leftSpread, double rightSpread);

  /**
   * Constructor.
   *
   * @param leftMean      left mean
   * @param rightMean     right mean
   * @param leftSpread    left spread
   * @param rightSpread   right spread
   */
  TriangularFuzzyInterval(double leftMean, double rightMean,
                          double leftSpread, double rightSpread);

  /**
   * Copy constructor.
   *
   * @param other   other triangular fuzzy interval
   */
  TriangularFuzzyInterval(const TriangularFuzzyInterval& other);

  /**
   * Destructor.
   */
  ~TriangularFuzzyInterval() override;

  /**
   * Evaluate the membership function.
   *
   * @param x   \f$x \in X\f$
   * @return    \f$\mu_{\tilde{x}}(x) \in [0, 1]\f$
   */
  double evaluateMembershipFunction(double x) const override;

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
   * @return  left mean
   */
  double getLeftMean() const;

  /**
   * @return  right mean
   */
  double getRightMean() const;

  /**
   * @return  left spread
   */
  double getLeftSpread() const;

  /**
   * @return  right spread
   */
  double getRightSpread() const;

 protected:
  /// left mean
  double leftMean;
  /// right mean
  double rightMean;
  /// left spread
  double leftSpread;
  /// right spread
  double rightSpread;
};

}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_FUZZY_TRIANGULARFUZZYINTERVAL_HPP */
