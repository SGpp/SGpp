// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_FUZZY_FUZZYEXTENSIONPRINCIPLE_HPP
#define SGPP_OPTIMIZATION_FUZZY_FUZZYEXTENSIONPRINCIPLE_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/base/function/scalar/ScalarFunction.hpp>
#include <sgpp/optimization/fuzzy/FuzzyInterval.hpp>

#include <vector>

namespace sgpp {
namespace optimization {

/**
 * Abstract class for Zadeh's fuzzy extension principle to propagate fuzzy input
 * uncertainties through a function to obtain a fuzzy output uncertainty.
 * Subclasses have to specify how the output confidence intervals are determined.
 *
 * Literature: W. Andreas Klimke. Uncertainty Modeling using Fuzzy Arithmetic and
 * Sparse Grids. PhD thesis, University of Stuttgart, IANS, 2006.
 */
class FuzzyExtensionPrinciple {
 public:
  /// default number of \f$\alpha\f$ segments
  static const size_t DEFAULT_NUMBER_OF_ALPHA_SEGMENTS = 10;

  /**
   * Constructor.
   *
   * @param f                       function through which to propagate the uncertainties
   * @param numberOfAlphaSegments   number of \f$\alpha\f$ segments
   */
  explicit FuzzyExtensionPrinciple(
      const base::ScalarFunction& f,
      size_t numberOfAlphaSegments = DEFAULT_NUMBER_OF_ALPHA_SEGMENTS);

  /**
   * Copy constructor.
   *
   * @param other   other fuzzy extension principle
   */
  FuzzyExtensionPrinciple(const FuzzyExtensionPrinciple& other);

  /**
   * Destructor.
   */
  virtual ~FuzzyExtensionPrinciple();

  /**
   * Apply the fuzzy extension principle to fuzzy input intervals.
   *
   * @param xFuzzy  vector of fuzzy input intervals
   * @return        fuzzy output interval; plain pointer due to the Python interface;
   *                wrap return value with smart pointer
   */
  FuzzyInterval* apply(const std::vector<const FuzzyInterval*>& xFuzzy);

  /**
   * @return  number of \f$\alpha\f$ segments
   */
  size_t getNumberOfAlphaSegments() const;

  /**
   * @param numberOfAlphaSegments   number of \f$\alpha\f$ segments
   */
  void setNumberOfAlphaSegments(size_t numberOfAlphaSegments);

  /**
   * @return  vector of \f$\alpha\f$ levels (after <tt>apply</tt> call)
   */
  const base::DataVector& getAlphaLevels() const;

  /**
   * @return  vector of lower bounds of input confidence intervals
   *          (after <tt>apply</tt> call)
   */
  const std::vector<base::DataVector>& getOptimizationDomainsLowerBounds() const;

  /**
   * @return  vector of upper bounds of input confidence intervals
   *          (after <tt>apply</tt> call)
   */
  const std::vector<base::DataVector>& getOptimizationDomainsUpperBounds() const;

  /**
   * @return  vector of minimum points (after <tt>apply</tt> call)
   */
  const std::vector<base::DataVector>& getMinimumPoints() const;

  /**
   * @return  vector of minimum function values (after <tt>apply</tt> call)
   */
  const base::DataVector& getMinimumValues() const;

  /**
   * @return  vector of maximum points (after <tt>apply</tt> call)
   */
  const std::vector<base::DataVector>& getMaximumPoints() const;

  /**
   * @return  vector of maximum function values (after <tt>apply</tt> call)
   */
  const base::DataVector& getMaximumValues() const;

  /**
   * Pure virtual method for cloning the fuzzy extension principle.
   * It should generate a pointer to the cloned object and
   * it's used for parallel computations
   * (the optimizeForSingleAlphaLevel() method might not be thread-safe).
   *
   * @param[out]  clone pointer to cloned object
   */
  virtual void clone(std::unique_ptr<FuzzyExtensionPrinciple>& clone) const = 0;

 protected:
  /// function through which to propagate the uncertainties
  std::unique_ptr<base::ScalarFunction> f;
  /// number of \f$\alpha\f$ segments
  size_t m;
  /// number of \f$\alpha\f$ segments
  base::DataVector alphaLevels;
  /// vector of lower bounds of input confidence intervals (after <tt>apply</tt> call)
  std::vector<base::DataVector> optimizationDomainsLowerBounds;
  /// vector of upper bounds of input confidence intervals (after <tt>apply</tt> call)
  std::vector<base::DataVector> optimizationDomainsUpperBounds;
  /// vector of minimum points (after <tt>apply</tt> call)
  std::vector<base::DataVector> minimumPoints;
  /// vector of minimum function values (after <tt>apply</tt> call)
  base::DataVector minimumValues;
  /// vector of maximum points (after <tt>apply</tt> call)
  std::vector<base::DataVector> maximumPoints;
  /// vector of maximum function values (after <tt>apply</tt> call)
  base::DataVector maximumValues;

  /**
   * Custom preparation method that is called before the parallelized
   * optimizeForSingleAlphaLevel calls. Here empty, but can be overridden by subclasses.
   */
  virtual void prepareApply();

  /**
   * Pure virtual method for solving the minimization/maximization problem
   * for a single \f$\alpha\f$ level.
   *
   * @param[in]   j             index of \f$\alpha\f$ level
   * @param[out]  minimumPoint  minimum point
   * @param[out]  minimumValue  minimum function value
   * @param[out]  maximumPoint  maximum point
   * @param[out]  maximumValue  maximum function value
   */
  virtual void optimizeForSingleAlphaLevel(
      size_t j, base::DataVector& minimumPoint, double& minimumValue,
      base::DataVector& maximumPoint, double& maximumValue) = 0;
};

}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_FUZZY_FUZZYEXTENSIONPRINCIPLE_HPP */
