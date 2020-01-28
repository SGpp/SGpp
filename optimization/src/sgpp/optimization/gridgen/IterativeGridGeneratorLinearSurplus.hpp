// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/gridgen/IterativeGridGenerator.hpp>

#include <cstddef>

namespace sgpp {
namespace optimization {

/**
 * Iterative grid generation based on linear surplusses.
 * In each iteration, the fraction of \f$gamma\f$
 * (e.g. \f$\gamma = 0.2\f$ means 20%)
 * of the grid points with the largest hierarchical linear
 * surplusses are refined.
 */
class IterativeGridGeneratorLinearSurplus : public IterativeGridGenerator {
 public:
  /// default adaptivity
  static constexpr double DEFAULT_ADAPTIVITY = 0.2;
  /// default level of initial regular sparse grid
  static const base::level_t DEFAULT_INITIAL_LEVEL = 3;

  /**
   * Constructor.
   * Do not destruct the grid before this object!
   *
   * @param f             objective function
   * @param grid          grid (should be empty)
   * @param N             maximal number of grid points
   * @param adaptivity    adaptivity between 0 and 1
   * @param initialLevel  level of initial regular sparse grid
   */
  IterativeGridGeneratorLinearSurplus(base::ScalarFunction& f, base::Grid& grid, size_t N,
                                      double adaptivity = DEFAULT_ADAPTIVITY,
                                      base::level_t initialLevel = DEFAULT_INITIAL_LEVEL);

  /**
   * Destructor.
   */
  ~IterativeGridGeneratorLinearSurplus() override;

  /**
   * Generate the grid.
   *
   * @return true on success, otherwise false
   */
  bool generate() override;

  /**
   * @return            adaptivity between 0 and 1
   */
  double getAdaptivity() const;

  /**
   * @param adaptivity  adaptivity between 0 and 1
   */
  void setAdaptivity(double adaptivity);

  /**
   * @return              level of initial regular sparse grid
   */
  base::level_t getInitialLevel() const;

  /**
   * @param initialLevel  level of initial regular sparse grid
   */
  void setInitialLevel(base::level_t initialLevel);

 protected:
  /// pointer to a linear grid (of the same type as the "real" grid)
  std::unique_ptr<base::Grid> linearGrid;
  /// adaptivity
  double gamma;
  /// level of initial regular sparse grid
  base::level_t initialLevel;
};
}  // namespace optimization
}  // namespace sgpp
