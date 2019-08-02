// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_FUZZY_FUZZYINTERVAL_HPP
#define SGPP_OPTIMIZATION_FUZZY_FUZZYINTERVAL_HPP

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace optimization {

/**
 * Abstract class for a fuzzy interval.
 *
 * A fuzzy set is the graph \f$\tilde{x} = \{(x, \mu_{\tilde{x}}(x)) \mid x \in X\}\f$
 * of some function \f$\mu_{\tilde{x}}\colon X \to [0, 1]\f$ (membership function),
 * where usually \f$X = [0, 1]\f$ (since we are working on sparse grids here).
 * A fuzzy interval is convex if
 * \f$\min(\mu_{\tilde{x}}(a), \mu_{\tilde{x}}(c)) \le \mu_{\tilde{x}}(b)\f$
 * for all \f$a, b, c \in X\f$ with \f$a \le b \le c\f$.
 * A fuzzy interval is normalized if \f$\max \mu_{\tilde{x}} = 1\f$.
 * A fuzzy interval is a convex and normalized fuzzy set with piecewise
 * continuous membership function.
 *
 * For some given \f$\alpha \in [0, 1]\f$, the confidence interval of level
 * \f$\alpha\f$ or the \f$\alpha\f$-cut is defined as
 * \f$(\tilde{x})_\alpha = \{x \in X \mid \mu_{\tilde{x}}(x) \ge \alpha\}\f$
 * for \f$\alpha > 0\f$ and
 * \f$(\tilde{x})_\alpha = \mathrm{supp}(\mu_{\tilde{x}})\f$
 * for \f$\alpha = 0\f$.
 * The confidence intervals of fuzzy intervals are always nested closed intervals,
 * i.e., \f$(\tilde{x})_\alpha = [a, b]\f$ for some \f$a \le b\f$ and
 * \f$(\tilde{x})_{\alpha_1} \supset (\tilde{x})_{\alpha_2}\f$ for
 * \f$\alpha_1 \le \alpha_2\f$.
 */
class FuzzyInterval {
 public:
  /// mode to determine norms of the fuzzy interval
  enum class NormMode {
    /// determine norm of membership function
    ViaMembershipFunction,
    /// determine norm of size of confidence interval
    ViaConfidenceInterval,
  };

  /// default number of samples to compute norms
  const size_t DEFAULT_NUMBER_OF_INTEGRAL_SAMPLES = 10000;

  /**
   * Constructor. Needs the support of the fuzzy interval (which is
   * always a closed interval, so it suffices to supply lower and upper bound).
   *
   * @param supportLowerBound   lower bound of the support
   * @param supportUpperBound   upper bound of the support
   */
  FuzzyInterval(double supportLowerBound, double supportUpperBound);

  /**
   * Destructor.
   */
  virtual ~FuzzyInterval();

  /**
   * Pure virtual method to evaluate the membership function.
   *
   * @param x   \f$x \in X\f$
   * @return    \f$\mu_{\tilde{x}}(x) \in [0, 1]\f$
   */
  virtual double evaluateMembershipFunction(double x) const = 0;

  /**
   * Pure virtual method to evaluate the lower bound of a confidence interval,
   * which is always a closed interval \f$(\tilde{x})_\alpha = [a, b]\f$.
   *
   * @param alpha   \f$\alpha \in [0, 1]\f$
   * @return        \f$a \in X\f$
   */
  virtual double evaluateConfidenceIntervalLowerBound(double alpha) const = 0;

  /**
   * Pure virtual method to evaluate the upper bound of a confidence interval,
   * which is always a closed interval \f$(\tilde{x})_\alpha = [a, b]\f$.
   *
   * @param alpha   \f$\alpha \in [0, 1]\f$
   * @return        \f$b \in X\f$
   */
  virtual double evaluateConfidenceIntervalUpperBound(double alpha) const = 0;

  /**
   * Compute L1 norm of fuzzy interval.
   *
   * @param normMode  mode with which to compute the norm
   * @return          L1 norm
   */
  double computeL1Norm(NormMode normMode = NormMode::ViaMembershipFunction) const;

  /**
   * Compute L2 norm of fuzzy interval.
   *
   * @param normMode  mode with which to compute the norm
   * @return          L2 norm
   */
  double computeL2Norm(NormMode normMode = NormMode::ViaMembershipFunction) const;

  /**
   * Compute Linf norm of fuzzy interval.
   *
   * @param normMode  mode with which to compute the norm
   * @return          Linf norm
   */
  double computeLinfNorm(NormMode normMode = NormMode::ViaMembershipFunction) const;

  /**
   * Compute absolute L1 error to other fuzzy interval.
   *
   * @param other     other fuzzy interval
   * @param normMode  mode with which to compute the norms
   * @return          absolute L1 error
   */
  double computeL1Error(
      const FuzzyInterval& other,
      NormMode normMode = NormMode::ViaMembershipFunction) const;

  /**
   * Compute absolute L2 error to other fuzzy interval.
   *
   * @param other     other fuzzy interval
   * @param normMode  mode with which to compute the norms
   * @return          absolute L2 error
   */
  double computeL2Error(
      const FuzzyInterval& other,
      NormMode normMode = NormMode::ViaMembershipFunction) const;

  /**
   * Compute absolute Linf error to other fuzzy interval.
   *
   * @param other     other fuzzy interval
   * @param normMode  mode with which to compute the norms
   * @return          absolute Linf error
   */
  double computeLinfError(
      const FuzzyInterval& other,
      NormMode normMode = NormMode::ViaMembershipFunction) const;

  /**
   * Compute relative L1 error to other fuzzy interval.
   *
   * @param other     other fuzzy interval
   * @param normMode  mode with which to compute the norms
   * @return          relative L1 error
   */
  double computeRelativeL1Error(
      const FuzzyInterval& other,
      NormMode normMode = NormMode::ViaMembershipFunction) const;

  /**
   * Compute relative L2 error to other fuzzy interval.
   *
   * @param other     other fuzzy interval
   * @param normMode  mode with which to compute the norms
   * @return          relative L2 error
   */
  double computeRelativeL2Error(
      const FuzzyInterval& other,
      NormMode normMode = NormMode::ViaMembershipFunction) const;

  /**
   * Compute relative Linf error to other fuzzy interval.
   *
   * @param other     other fuzzy interval
   * @param normMode  mode with which to compute the norms
   * @return          relative Linf error
   */
  double computeRelativeLinfError(
      const FuzzyInterval& other,
      NormMode normMode = NormMode::ViaMembershipFunction) const;

  /**
   * @return lower bound of the support
   */
  double getSupportLowerBound() const;

  /**
   * @return upper bound of the support
   */
  double getSupportUpperBound() const;

  /**
   * @return number of samples to compute norms
   */
  size_t getNumberOfIntegralSamples() const;

  /**
   * @param numberOfIntegralSamples   number of samples to compute norms
   */
  void setNumberOfIntegralSamples(size_t numberOfIntegralSamples);

 protected:
  /// lower bound of the support
  double supportLowerBound;
  /// upper bound of the support
  double supportUpperBound;
  /// number of samples to compute norms
  size_t numberOfIntegralSamples;
};

}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_FUZZY_FUZZYINTERVAL_HPP */
