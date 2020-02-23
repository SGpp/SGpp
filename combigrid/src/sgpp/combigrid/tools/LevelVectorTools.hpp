// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/combigrid/LevelIndexTypes.hpp>
#include <sgpp/globaldef.hpp>

#include <vector>

namespace sgpp {
namespace combigrid {

/**
 * Various tools for LevelVectors, mainly for generating specific vectors of LevelVectors.
 */
class LevelVectorTools {
 public:
  LevelVectorTools() = delete;

  /**
   * @brief Get a hypercube of all level vectors between and including minLevel and maxLevel.
   *
   * @param maxLevel  the maximum level vector
   * @param minLevel  the minimum level vector, must not be larger than maxLevel in any dimension
   * @return a lexicographically ordered vector of level vectors in the hypercube
   */
  static std::vector<LevelVector> generateHyperCube(const LevelVector& maxLevel,
                                                    const LevelVector& minLevel);

  /**
   * @brief Enumerate all levels \f$\vec{\ell} \in \mathbb{N}_{\ge 0}^{dim}\f$
   * with \f$\sum_{d=1}^{dim} \ell_d = \mathrm{levelSum}\f$ and
   * \f$\ell_d \ge \textrm{minLevel}_d\f$ for all \f$d\f$.
   *
   * @param minLevel  the minimum level vector
   * @param levelSum  level sum
   * @return vector of levels of desired level sum
   */
  static std::vector<LevelVector> generateDiagonal(const LevelVector& minLevel, level_t levelSum);

  /**
   * @brief Enumerate all levels \f$\vec{\ell} \in \mathbb{N}_{\ge 0}^{dim}\f$
   * with \f$\sum_{d=1}^{dim} \ell_d = n\f$.
   *
   * @param dim       dimensionality
   * @param levelSum  level sum
   * @return vector of levels of desired level sum (with boundary)
   */
  static std::vector<LevelVector> generateDiagonalWithBoundary(size_t dim, level_t levelSum);

  /**
   * @brief Enumerate all levels \f$\vec{\ell} \in \mathbb{N}_{\ge 1}^{dim}\f$
   * with \f$\sum_{d=1}^{dim} \ell_d = n\f$.
   *
   * @param dim       dimensionality
   * @param levelSum  level sum
   * @return vector of levels of desired level sum (without boundary)
   */
  static std::vector<LevelVector> generateDiagonalWithoutBoundary(size_t dim, level_t levelSum);

  /**
   * @brief Make any level set a downward closed one.
   *
   * @param subspaceLevels      an arbitrary level set
   * @param lowestLevelVector   the minimum level vector, a vector of zeros or ones would be a
   *                                standard choice
   * @return the updated, downward closed level set containing subspaceLevels
   */
  static std::vector<LevelVector> makeDownwardClosed(const std::vector<LevelVector>& subspaceLevels,
                                                     LevelVector lowestLevelVector);

 protected:
  static std::vector<LevelVector> generateHyperCubeRecursiveLastDim(const LevelVector& maxLevel,
                                                                    const LevelVector& minLevel,
                                                                    const LevelVector& prefix);
  static std::vector<LevelVector> generateHyperCubeRecursive(const LevelVector& maxLevel,
                                                             const LevelVector& minLevel,
                                                             const LevelVector& prefix);

  static std::vector<LevelVector> generateDiagonalRecursive(const LevelVector& minLevel,
                                                            level_t minLevelSum,
                                                            level_t levelSum,
                                                            const LevelVector& suffix);
};

}  // namespace combigrid
}  // namespace sgpp
