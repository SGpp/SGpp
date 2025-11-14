// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/gridgen/IterativeGridGenerator.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

namespace sgpp {
namespace optimization {
/**
 * @brief Helper classes and functions for the IterativeGridGeneratorFullAdaptiveRitterNovak.
 */
namespace IGGFARNHelper {

/**
 * @brief Stores and normalizes the hyperparameters for the refinement criterion.
 */
class Hyperparameters {
 private:
  double levelDegree;
  double threshold;
  double extendedThreshold;
  double gradient;
  double secondGradient;
  double gaussianProcess;

 public:
  /**
   * @brief Constructor that takes individual weights for each criterion.
   * The weights are normalized to sum to 1.
   *
   * @param levelDegree       Weight for the level-based criterion.
   * @param threshold         Weight for the threshold criterion (value at parent).
   * @param extendedThreshold Weight for the extended threshold criterion (value at children).
   * @param gradient          Weight for the gradient-based criterion.
   * @param secondGradient    Weight for the second-order gradient criterion.
   * @param gaussianProcess   Weight for the Gaussian Process uncertainty criterion.
   */
  Hyperparameters(double levelDegree, double threshold, double extendedThreshold, double gradient,
                  double secondGradient, double gaussianProcess);

  /**
   * @brief Default constructor, initializes all weights to zero.
   */
  Hyperparameters();

  double getLevelDegree() const;
  double getThreshold() const;
  double getExtendedThreshold() const;
  double getGradient() const;
  double getSecondGradient() const;
  double getGaussianProcess() const;
};

/// Stream operator for printing hyperparameters.
std::ostream& operator<<(std::ostream& os, const Hyperparameters& hyperparameters);

/**
 * @brief A class to represent a value with uncertainty bounds, typically from sampling.
 */
class UncertaintyBound {
 private:
  double _lower;
  double _value;
  double _upper;
  double _coV;  // Coefficient of Variation

 public:
  /**
   * @brief Construct from explicit lower, value, upper, and CoV.
   *
   * @param lower Lower bound.
   * @param value The estimated value.
   * @param upper Upper bound.
   * @param coV   Coefficient of variation.
   */
  UncertaintyBound(double lower, double value, double upper, double coV);

  /**
   * @brief Default constructor, initializes to zero.
   */
  UncertaintyBound();

  /**
   * @brief Construct from a number of positive samples out of a total number of samples.
   * Calculates the value and confidence interval.
   *
   * @param samples Number of positive samples.
   * @param N       Total number of samples.
   */
  UncertaintyBound(const size_t samples, const size_t N);

  /**
   * @brief Divides this uncertainty bound by another. Used for relative error estimation.
   *
   * @param other The divisor.
   * @return A new UncertaintyBound representing the quotient.
   */
  UncertaintyBound div(const UncertaintyBound& other) const;

  double lower() const;
  double value() const;
  double upper() const;
  double coV() const;
};

/// Stream operator for printing uncertainty bounds.
std::ostream& operator<<(std::ostream& os, const UncertaintyBound& bounds);

/**
 * @brief Implements a stopping criterion based on the convergence of the surrogate model.
 *
 * It compares the current grid's surrogate with surrogates from previous iterations. If the
 * relative mismatch probability between them is consistently below a tolerance, it signals to stop.
 */
class StoppingCriterion {
 protected:
  static const size_t stopping_criterion_count_difference = 50;
  static constexpr double stopping_criterion_tolerance = 0.025;

  const std::vector<sgpp::base::DataVector>& monteCarloPoints;
  std::vector<
      std::tuple<size_t, std::vector<double>, std::unique_ptr<base::Grid>, base::DataVector>>
      last_grid_values;

  /**
   * @brief Measures the failure and mismatch probabilities between two functions.
   *
   * @param points             Points to sample from.
   * @param surrogate_values   Cache for values of the first function.
   * @param original_values    Cache for values of the second function.
   * @param surrogate_f        First function (e.g., current surrogate).
   * @param original_f         Second function (e.g., previous surrogate).
   * @param early_stop_criterion A predicate to decide if sampling can stop early.
   * @return A tuple containing (number of samples used, failure probability, mismatch
   * probability).
   */
  std::tuple<size_t, UncertaintyBound, UncertaintyBound> measure_error(
      const std::vector<base::DataVector>& points, std::vector<double>& surrogate_values,
      std::vector<double>& original_values,
      const std::function<double(const sgpp::base::DataVector&)>& surrogate_f,
      const std::function<double(const sgpp::base::DataVector&)>& original_f,
      const std::function<bool(const UncertaintyBound&, const UncertaintyBound&, size_t)>&
          early_stop_criterion);

 public:
  /**
   * @brief Constructor.
   *
   * @param monteCarloPoints Points to use for validation.
   */
  explicit StoppingCriterion(const std::vector<sgpp::base::DataVector>& monteCarloPoints);

  /**
   * @brief Destructor.
   */
  ~StoppingCriterion();

  /**
   * @brief Checks if the grid generation should stop.
   *
   * @param grid  The current grid.
   * @param alpha The current surplus vector.
   * @return True if the stopping criterion is met, false otherwise.
   */
  bool shouldStop(base::Grid& grid, const base::DataVector& alpha);
};
}  // namespace IGGFARNHelper

/**
 * @brief Iterative grid generation based on a modified, fully adaptive Ritter/Novak's refinement
 * criterion.
 *
 * This generator uses a combination of weighted criteria to decide which grid point to refine.
 * The criteria include the grid point's level, function value (threshold), gradient, and more.
 * A large number of hyperparameters allow these criteria to be weighted differently.
 *
 * @warning This class uses HashRefinementMultiple, so it generates grids that may not meet the
 * "hierarchical ancestors" requirement.
 *
 * Literature: Julius von Smercek, Bachelor Thesis ``Critical Event Estimation using Gradient-based
 * Refinement Methods for Adaptive Sparse Grids'', 2025
 *
 * @see HashRefinementMultiple
 */
class IterativeGridGeneratorFullAdaptiveRitterNovak : public IterativeGridGenerator {
 public:
  /// Default level of the initial regular sparse grid.
  static const base::level_t DEFAULT_INITIAL_LEVEL = 3;
  /// Default maximal level of grid points.
  static const base::level_t DEFAULT_MAX_LEVEL = 20;
#ifdef USE_LIBGP
  static const bool LIBGP_AVAILABLE = true;
#else
  static const bool LIBGP_AVAILABLE = false;
#endif

  /**
   * @brief Constructor.
   * Do not destruct the grid object before this generator object!
   *
   * @param f                 Objective function.
   * @param grid              Grid object (should be empty).
   * @param N                 Maximal number of grid points.
   * @param hyperparameters   Hyperparameters for the adaptive weight function.
   * @param domain            Domain for each dimension
   * @param refineLeftOrRight Whether to only refine the left or right child in each dimension,
   * activating the "full adaptive" approach.
   * @param initialLevel      Level of the initial regular sparse grid.
   * @param maxLevel          Maximal level permitted for any grid point.
   */
  IterativeGridGeneratorFullAdaptiveRitterNovak(
      base::ScalarFunction& f, base::Grid& grid, size_t N,
      const IGGFARNHelper::Hyperparameters& hyperparameters,
      const std::vector<std::pair<double, double>>& domain, bool refineLeftOrRight = true,
      const base::level_t initialLevel = DEFAULT_INITIAL_LEVEL,
      const base::level_t maxLevel = DEFAULT_MAX_LEVEL);

  /**
   * @brief Destructor.
   */
  ~IterativeGridGeneratorFullAdaptiveRitterNovak() override;

  /**
   * @brief Generates the grid adaptively until the maximum number of points `N` is reached
   * or a stopping criterion is met.
   *
   * @return True on success, false otherwise.
   */
  bool generate() override;

  /**
   * @return The currently used hyperparameters.
   */
  const IGGFARNHelper::Hyperparameters& getHyperparameters() const;

  /**
   * @param hyperparameters New hyperparameters for the adaptive weight function.
   */
  void setHyperparameters(const IGGFARNHelper::Hyperparameters& hyperparameters);

  /**
   * @return Level of the initial regular sparse grid.
   */
  base::level_t getInitialLevel() const;

  /**
   * @param initialLevel New level for the initial regular sparse grid.
   */
  void setInitialLevel(base::level_t initialLevel);

  /**
   * @return Maximal permitted level for grid points.
   */
  base::level_t getMaxLevel() const;

  /**
   * @param maxLevel New maximal level for grid points.
   */
  void setMaxLevel(base::level_t maxLevel);

  /**
   * @brief Activates the stopping criterion.
   * The generation will stop if the surrogate model's predictions stabilize between iterations.
   *
   * @param validationPoints A set of points used to compare surrogate models.
   */
  void activateStoppingCriterion(const std::vector<sgpp::base::DataVector>& validationPoints);

  /**
   * @brief Disables the stopping criterion.
   */
  void disableStoppingCriterion();

  /**
   * @return True if the stopping criterion is enabled, false otherwise.
   */
  bool isStoppingCriterionEnabled() const;

  /**
   * @return The coefficients (surpluses) of the final sparse grid interpolant.
   */
  const base::DataVector& getAlpha() const;

 protected:
  /// Hyperparameters for the adaptive weight function.
  IGGFARNHelper::Hyperparameters hyperparameters;
  /// Domain bounds for each dimension.
  std::vector<std::pair<double, double>> domain;
  /// Level of the initial regular sparse grid.
  base::level_t initialLevel;
  /// Maximal level of grid points.
  base::level_t maxLevel;
  /// Validation points for the stopping criterion. Empty if disabled.
  std::vector<sgpp::base::DataVector> stoppingCriterionValidationPoints;
  /// Whether to refine only the left or right child in each dimension.
  bool refineLeftOrRight;
  /// Coefficients of the sparse grid surrogate.
  base::DataVector alpha;
};
}  // namespace optimization
}  // namespace sgpp
