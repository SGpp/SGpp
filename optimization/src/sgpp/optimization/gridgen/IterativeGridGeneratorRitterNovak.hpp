// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/gridgen/IterativeGridGenerator.hpp>

#include <cstddef>

namespace sgpp {
namespace optimization {

/**
 * Iterative grid generation based on Ritter/Novak's refinement criterion.
 * Caution: This class uses HashRefinementMultiple, so it generates grids
 * that don't meet the "hierarchical ancestors" requirement!
 *
 * Literature: Erich Novak, Klaus Ritter. Global Optimization Using
 * Hyperbolic Cross Points.
 * In: Christodoulos A. Floudas, Panos M. Pardalos (eds.): State of the
 * Art in Global Optimization,
 * Computational Methods and Applications, Vol. 7. Springer 1996.
 * DOI: 10.1007/978-1-4613-3437-8_2
 *
 * @see HashRefinementMultiple
 */
class IterativeGridGeneratorRitterNovak : public IterativeGridGenerator {
 public:
  /// default adaptivity
  static constexpr double DEFAULT_ADAPTIVITY = 0.85;
  /// default level of initial regular sparse grid
  static const base::level_t DEFAULT_INITIAL_LEVEL = 3;
  /// default maximal level of grid points
  static const base::level_t DEFAULT_MAX_LEVEL = 20;

  /// exponentiation methods
  enum PowMethod { STD_POW, FAST_POW };

  /**
   * Constructor.
   * Do not destruct the grid before this object!
   *
   * @param f             objective function
   * @param grid          grid (should be empty)
   * @param N             maximal number of grid points
   * @param adaptivity    adaptivity between 0 and 1
   * @param initialLevel  level of initial regular sparse grid
   * @param maxLevel      maximal level of grid points
   * @param powMethod     exponentiation method
   *                      (fastPow is faster than std::pow,
   *                      but only approximative)
   */
  IterativeGridGeneratorRitterNovak(base::ScalarFunction& f, base::Grid& grid, size_t N,
                                    double adaptivity = DEFAULT_ADAPTIVITY,
                                    base::level_t initialLevel = DEFAULT_INITIAL_LEVEL,
                                    base::level_t maxLevel = DEFAULT_MAX_LEVEL,
                                    PowMethod powMethod = STD_POW);

  /**
   * Destructor.
   */
  ~IterativeGridGeneratorRitterNovak() override;

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

  /**
   * @return          maximal level of grid points
   */
  base::level_t getMaxLevel() const;

  /**
   * @param maxLevel  maximal level of grid points
   */
  void setMaxLevel(base::level_t maxLevel);

  /**
   * @return          exponentiation method
   */
  PowMethod getPowMethod() const;

  /**
   * @param powMethod exponentiation method
   */
  void setPowMethod(PowMethod powMethod);

 protected:
  /// adaptivity
  double gamma;
  /// level of initial regular sparse grid
  base::level_t initialLevel;
  /// maximal level of grid points
  base::level_t maxLevel;
  /// exponentiation method
  PowMethod powMethod;
};
}  // namespace optimization
}  // namespace sgpp
