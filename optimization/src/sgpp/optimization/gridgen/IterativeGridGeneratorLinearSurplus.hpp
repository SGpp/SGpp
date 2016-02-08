// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_GRIDGEN_ITERATIVEGRIDGENERATORLINEARSURPLUS_HPP
#define SGPP_OPTIMIZATION_GRIDGEN_ITERATIVEGRIDGENERATORLINEARSURPLUS_HPP

#include <cstddef>

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/gridgen/IterativeGridGenerator.hpp>

namespace SGPP {
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
  static constexpr float_t DEFAULT_ADAPTIVITY = 0.2;
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
  IterativeGridGeneratorLinearSurplus(
    ScalarFunction& f,
    base::Grid& grid,
    size_t N,
    float_t adaptivity = DEFAULT_ADAPTIVITY,
    base::level_t initialLevel = DEFAULT_INITIAL_LEVEL);

  /**
   * Destructor.
   */
  virtual ~IterativeGridGeneratorLinearSurplus() override;

  /**
   * Generate the grid.
   *
   * @return true on success, otherwise false
   */
  virtual bool generate() override;

  /**
   * @return            adaptivity between 0 and 1
   */
  float_t getAdaptivity() const;

  /**
   * @param adaptivity  adaptivity between 0 and 1
   */
  void setAdaptivity(float_t adaptivity);

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
  float_t gamma;
  /// level of initial regular sparse grid
  base::level_t initialLevel;
};

}
}

#endif /* SGPP_OPTIMIZATION_GRIDGEN_ITERATIVEGRIDGENERATORLINEARSURPLUS_HPP */
