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

/**
 * Quasi-Gaussian fuzzy number. A fuzzy number is a fuzzy interval where
 * \f$\{x \in X \mid \mu_{\tilde{x}}(x) = 1\} = \{a\}\f$ for some \f$a \in X\f$.
 * Quasi-Gaussian fuzzy numbers have a cut-off Gaussian function as membership
 * function, which is parametrized by its mean, the standard deviation, and the
 * cut-off point.
 */
class QuasiGaussianFuzzyNumber : public FuzzyIntervalViaMembershipFunction {
 public:
  /**
   * Constructor.
   *
   * @param mean    mean
   * @param stdev   standard deviation
   * @param cutoff  cut-off point
   */
  QuasiGaussianFuzzyNumber(double mean, double stdev, double cutoff);

  /**
   * Copy constructor.
   *
   * @param other   other quasi-Gaussian fuzzy number
   */
  QuasiGaussianFuzzyNumber(const QuasiGaussianFuzzyNumber& other);

  /**
   * Destructor.
   */
  ~QuasiGaussianFuzzyNumber() override;

  /**
   * Evaluate the membership function.
   *
   * @param x   \f$x \in X\f$
   * @return    \f$\mu_{\tilde{x}}(x) \in [0, 1]\f$
   */
  double evaluateMembershipFunction(double x) const override;

  /**
   * @return  mean
   */
  double getMean() const;

  /**
   * @return  standard deviation
   */
  double getStdev() const;

  /**
   * @return  cut-off point
   */
  double getCutoff() const;

 protected:
  /// mean
  double mean;
  /// standard deviation
  double stdev;
  /// cut-off point
  double cutoff;
};

}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_FUZZY_QUASIGAUSSIANFUZZYNUMBER_HPP */
