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

/**
 * Fuzzy interval by piecewise linear interpolation of sample points
 * \f$(x_i, \alpha_i)\f$ of the membership function \f$\mu_{\tilde{x}}\f$, i.e.,
 * \f$\mu_{\tilde{x}}(x_i) = \alpha_i\f$ (\f$i = 1, \dotsc, n\f$).
 *
 * The data must fulfill the following:
 * - \f$x_i < x_{i+1}\f$ for \f$i = 1, \dotsc, n-1\f$
 * - \f$\alpha_1 = 0 = \alpha_n\f$
 * - There are some \f$1 < j \le k < n\f$ such that:
 *   - \f$\alpha_i = 1\f$ for \f$i = j, \dotsc, k\f$
 *   - \f$\alpha_i < \alpha_{i+1}\f$ for \f$i = 1, \dotsc, j-1\f$
 *   - \f$\alpha_i > \alpha_{i+1}\f$ for \f$i = k, \dotsc, n-1\f$
 */
class InterpolatedFuzzyInterval : public FuzzyIntervalViaMembershipFunction {
 public:
  /**
   * Compute the lower bound of the core
   * (i.e., the area where the membership function equals 1).
   *
   * @param xData       \f$x\f$ data of sample points
   * @param alphaData   \f$\alpha\f$ data of sample points
   */
  static double getCoreLowerBound(const base::DataVector& xData,
                                  const base::DataVector& alphaData);

  /**
   * Compute the upper bound of the core
   * (i.e., the area where the membership function equals 1).
   *
   * @param xData       \f$x\f$ data of sample points
   * @param alphaData   \f$\alpha\f$ data of sample points
   */
  static double getCoreUpperBound(const base::DataVector& xData,
                                  const base::DataVector& alphaData);

  /**
   * Try to cast a FuzzyInterval to an InterpolatedFuzzyInterval
   * (<tt>dynamic_cast</tt>), needed for the Python interface.
   *
   * @param fuzzyInterval   FuzzyInterval to cast
   * @return                pointer to InterpolatedFuzzyInterval if succesful,
   *                        nullptr otherwise
   */
  static InterpolatedFuzzyInterval* tryDowncast(FuzzyInterval& fuzzyInterval);

  /**
   * Constructor.
   *
   * @param xData       \f$x\f$ data of sample points
   * @param alphaData   \f$\alpha\f$ data of sample points
   */
  InterpolatedFuzzyInterval(const base::DataVector& xData, const base::DataVector& alphaData);

  /**
   * Copy constructor.
   *
   * @param other   other interpolated fuzzy interval
   */
  InterpolatedFuzzyInterval(const InterpolatedFuzzyInterval& other);

  /**
   * Destructor.
   */
  ~InterpolatedFuzzyInterval() override;

  /**
   * Evaluate the membership function.
   *
   * @param x   \f$x \in X\f$
   * @return    \f$\mu_{\tilde{x}}(x) \in [0, 1]\f$
   */
  double evaluateMembershipFunction(double x) const override;

  /**
   * @return  \f$x\f$ data of sample points
   */
  const base::DataVector& getXData() const;

  /**
   * @return  \f$\alpha\f$ data of sample points
   */
  const base::DataVector& getAlphaData() const;

 protected:
  /// \f$x\f$ data of sample points
  base::DataVector xData;
  /// \f$\alpha\f$ data of sample points
  base::DataVector alphaData;
};

}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_FUZZY_INTERPOLATEDFUZZYINTERVAL_HPP */

