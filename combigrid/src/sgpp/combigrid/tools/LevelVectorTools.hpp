// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/combigrid/LevelIndexTypes.hpp>
#include <sgpp/globaldef.hpp>

#include <algorithm>
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
   * @brief Enumerate all levels in the hypercube between and including minLevel and maxLevel.
   *
   * @param minLevel  the minimum level vector
   * @param maxLevel  the maximum level vector
   * @return a colexicographically ordered vector of level vectors in the hypercube
   */
  static std::vector<LevelVector> generateHyperCube(const LevelVector& minLevel,
                                                    const LevelVector& maxLevel);

  /**
   * @brief Enumerate all levels in the hypercube between and including \f$(0, \dotsc, 0)\f$ and
   * maxLevel.
   *
   * @param maxLevel  the maximum level vector
   * @return a colexicographically ordered vector of level vectors in the hypercube
   */
  static std::vector<LevelVector> generateHyperCubeWithBoundary(const LevelVector& maxLevel);

  /**
   * @brief Enumerate all levels in the hypercube between and including \f$(1, \dotsc, 1)\f$ and
   * maxLevel.
   *
   * @param maxLevel  the maximum level vector
   * @return a colexicographically ordered vector of level vectors in the hypercube
   */
  static std::vector<LevelVector> generateHyperCubeWithoutBoundary(const LevelVector& maxLevel);

  /**
   * @brief Enumerate all levels \f$\vec{\ell} \in \mathbb{N}_{\ge 0}^{dim}\f$
   * with \f$\sum_{d=1}^{dim} \ell_d = \mathrm{levelSum}\f$ and
   * \f$\ell_d \ge \textrm{minLevel}_d\f$ for all \f$d\f$.
   *
   * @param minLevel  the minimum level vector
   * @param levelSum  level sum
   * @return a colexicographically ordered vector of levels of desired level sum
   */
  static std::vector<LevelVector> generateDiagonal(const LevelVector& minLevel, level_t levelSum);

  /**
   * @brief Enumerate all levels \f$\vec{\ell} \in \mathbb{N}_{\ge 0}^{dim}\f$
   * with \f$\sum_{d=1}^{dim} \ell_d = n\f$.
   *
   * @param dim       dimensionality
   * @param levelSum  level sum
   * @return a colexicographically ordered vector of levels of desired level sum (with boundary)
   */
  static std::vector<LevelVector> generateDiagonalWithBoundary(size_t dim, level_t levelSum);

  /**
   * @brief Enumerate all levels \f$\vec{\ell} \in \mathbb{N}_{\ge 1}^{dim}\f$
   * with \f$\sum_{d=1}^{dim} \ell_d = n\f$.
   *
   * @param dim       dimensionality
   * @param levelSum  level sum
   * @return a colexicographically ordered vector of levels of desired level sum (without boundary)
   */
  static std::vector<LevelVector> generateDiagonalWithoutBoundary(size_t dim, level_t levelSum);

  /**
   * @brief Make any level set a downward closed one.
   *
   * @param minLevel        the minimum level vector, a vector of zeros or ones would be a
   *                            standard choice
   * @param subspaceLevels  an arbitrary level set
   * @return the updated, downward closed level set containing subspaceLevels
   */
  static std::vector<LevelVector> makeDownwardClosed(
      LevelVector minLevel, const std::vector<LevelVector>& subspaceLevels);

  static void sort(std::vector<LevelVector>& levels);

  class Hash {
   public:
    size_t operator()(const LevelVector& level) const {
      std::hash<level_t> hasher;
      size_t seed = 0;
      for (level_t l : level) {
        seed ^= hasher(l) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      }
      return seed;
    }
  };

  static bool compareLower(const LevelVector& a, const LevelVector& b) {
    for (size_t dim = 0; dim < std::min(a.size(), b.size()); ++dim) {
      const size_t i = a.size() - dim - 1;
      const size_t j = b.size() - dim - 1;

      if (a[i] < b[j]) {
        return true;
      } else if (a[i] > b[j]) {
        return false;
      }
    }

    return (a.size() < b.size());
  }

 protected:
  /**
   * @brief Enumerate all levels in the hypercube between and including minLevel and maxLevel
   * (internal recursive function).
   *
   * @param[in] minLevel  the minimum level vector
   * @param[in] maxLevel  the maximum level vector
   * @param[in] curLevel  the current level vector, \c curLevel[d] won't be modified anymore for
   *                      \f$d >= \mathrm{curDim}\f$
   * @param[in] curDim    the dimension of the remaining hypercube
   * @param[out] result   a colexicographically ordered vector of level vectors in the hypercube
   */
  static void generateHyperCubeRecursive(const LevelVector& minLevel, const LevelVector& maxLevel,
                                         LevelVector& curLevel, size_t curDim,
                                         std::vector<LevelVector>& result);

  /**
   * @brief Enumerate all levels \f$\vec{\ell} \in \mathbb{N}_{\ge 0}^{dim}\f$
   * with \f$\sum_{d=1}^{dim} \ell_d = \mathrm{levelSum}\f$ and
   * \f$\ell_d \ge \textrm{minLevel}_d\f$ for all \f$d\f$ (internal recursive function).
   *
   * @param[in] minLevel      the minimum level vector
   * @param[in] minLevelSum   level sum of \c minLevel
   * @param[in] levelSum      level sum
   * @param[in] curLevel      the current level vector, \c curLevel[d] won't be modified anymore for
   *                              \f$d >= \mathrm{curDim}\f$
   * @param[in] curDim        the dimension of the remaining diagonal
   * @param[out] result       a colexicographically ordered vector of levels of desired level sum
   */
  static void generateDiagonalRecursive(const LevelVector& minLevel, level_t minLevelSum,
                                        level_t levelSum, LevelVector& curLevel, size_t curDim,
                                        std::vector<LevelVector>& result);
};

}  // namespace combigrid
}  // namespace sgpp
