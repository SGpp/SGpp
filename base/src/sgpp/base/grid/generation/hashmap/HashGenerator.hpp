// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef HASHGENERATOR_HPP
#define HASHGENERATOR_HPP

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/exception/generation_exception.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/globaldef.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iterator>
#include <set>
#include <unordered_set>
#include <vector>

namespace sgpp {
namespace base {

/**
 * This class provides the generation functionality of sparse grids
 * based on hashmaps.
 *
 * Grids with and without boundaries are supported.
 *
 * For boundary grids two cases are supported:
 *
 * 1. the classic sparse grid with level 0 and a diagonal
 * cut through the sub space scheme.
 *
 * 2. a modified boundary grid with level 0 and a pentagon cut
 * trough the sub space scheme.
 *
 * Furthermore, the creation of full grids (in the hierarchical basis) is supported.
 */
class HashGenerator {
 public:
  /**
   * Generates a regular sparse grid of level levels, without boundaries
   *
   * @param storage Hashmap that stores the grid points
   * @param level Grid level (non-negative value)
   * @param T modifier for subgrid selection, T = 0 implies standard sparse grid.
   *        For further information see Griebel and Knapek's paper
   *        optimized tensor-product approximation spaces.
   *        The effect of T can be seen in: \image html generalisedGrid.svg
   */
  void regular(GridStorage& storage, level_t level, double T = 0) {
    if (storage.getSize() > 0) {
      throw generation_exception("storage not empty");
    }
    this->regular_iter(storage, level, T);
  }

  /**
   * @brief Generates a regular sparse grid of level level, without boundaries.
   * The resulting grid only contains interactions that are in the vector terms.
   * @param storage Hashmap that stores the grid points
   * @param level Grid level (non-negative value)
   * @param terms controls the desired interaction terms.
   * For example, if we want to include grid points that model an
   * interaction between the first and the second predictor, we would
   * include the vector [1,2] in terms.
   * @param T modifier for subgrid selection, T = 0 implies standard sparse grid.
   *        For further information see Griebel and Knapek's paper
   *        optimized tensor-product approximation spaces.
   */
  void regularInter(GridStorage& storage, level_t level,
                    const std::set<std::set<size_t>>& terms, double T = 0) {
    if (storage.getSize() > 0) {
      throw generation_exception("storage not empty");
    }
    auto interset = std::unordered_set<std::vector<bool>>();
    const auto dim = storage.getDimension();
    for (const auto& interaction : terms) {
      auto term = std::vector<bool>(dim, false);
      for (const auto i : interaction) {
        term[i] = true;
      }
      interset.insert(term);
    }

    this->regular_inter(storage, level, interset, T);
  }

  void regular_inter(GridStorage& storage, level_t level,
                     const std::unordered_set<std::vector<bool>>& terms, double T = 0) {
    if (storage.getSize() > 0) {
      throw generation_exception("storage not empty");
    }
    this->regular_inter_iter(storage, level, terms, T);
  }

  /**
  * Generates a regular sparse grid of level levels, without boundaries
  * where dimensions are splitted into a groups with only certain number
  * of dimensions completely connected in a clique
  *
  * @param storage Hashmap that stores the grid points
  * @param level Grid level (non-negative value)
  * @param clique_size number of dimensions in a clique
  * @param T modifier for subgrid selection, T = 0 implies standard sparse grid.
  *        For further information see Griebel and Knapek's paper
  *        optimized tensor-product approximation spaces
  */
  void cliques(GridStorage& storage, level_t level, size_t clique_size, double T = 0) {
    if (storage.getSize() > 0) {
      throw generation_exception("storage not empty");
    }

    if (storage.getDimension() < clique_size) {
      throw generation_exception("clique size should be not greater than grid dimension");
    }

    this->cliques_iter(storage, level, clique_size, T);
  }

  /**
   * Generates a full grid of level @p level, without boundaries.
   *
   * @param storage Hashmap that stores the grid points
   * @param level Grid level (non-negative value)
   */
  void full(GridStorage& storage, level_t level) {
    if (storage.getSize() > 0) {
      throw generation_exception("storage not empty");
    }

    this->createFullGridIterative(storage, level);
  }
  /**
   * Generates an anisotropic full grid of the level vector
   *
   * @param storage Hashmap that stores the grid points
   * @param dimlevels grid level vector (vector of non-negative values)
   */
  void anisotropicFull(GridStorage& storage, std::vector<size_t>& dimlevels) {
    if (storage.getSize() > 0) {
      throw generation_exception("storage not empty");
    }

    this->createAnisotropicFullGrid(storage, dimlevels);
  }

  /**
   * Generates a full grid of level @p level, with boundary grid points.
   *
   * @param storage Hashmap that stores the grid points
   * @param level Grid level (non-negative value)
   */
  void fullWithBoundary(GridStorage& storage, level_t level) {
    if (storage.getSize() > 0) {
      throw generation_exception("storage not empty");
    }

    this->createFullGridTruncatedIterative(storage, level);
  }

  /**
   * Generates a regular sparse grid of level levels with boundaries
   *
   * @param storage Hashmap, that stores the grid points
   * @param level maximum level of the sparse grid (non-negative value)
   * @param boundaryLevel level at which the boundary points should be
   *                      inserted
   */
  void regularWithBoundaries(GridStorage& storage, level_t level, level_t boundaryLevel = 1) {
    if (storage.getSize() > 0) {
      throw generation_exception("storage not empty");
    }

    GridPoint point(storage.getDimension());

    if (boundaryLevel >= 1) {
      this->regular_boundary_truncated_iter(storage, level, boundaryLevel);
    } else {
      /* new grid generation
       *
       * for all level the same calculation of the level sum is implemented:
       * |l| <= n
       */
      for (size_t d = 0; d < storage.getDimension(); d++) {
        point.push(d, 0, 0, false);
      }

      this->boundaries_rec(storage, point, storage.getDimension() - 1, 0, level);
    }
  }

  /**
   * Generates a regular sparse grid of level levels with boundaries
   *
   * @param storage Hashmap, that stores the grid points
   * @param level maximum level of the sparse grid (non-negative value)
   * @param T modifier for subgrid selection, T = 0 implies standard sparse grid.
   *        For further information see Griebel and Knapek's paper
   *        optimized tensor-product approximation spaces
   */
  void regularWithPeriodicBoundaries(GridStorage& storage, level_t level, double T = 0) {
    if (storage.getSize() > 0) {
      throw generation_exception("storage not empty");
    }

    GridPoint point(storage.getDimension());

    if (level == 0) {
      throw generation_exception("Grid generation with maximum level 0 is not supported");
    } else {
      this->regular_periodic_boundary_iter(storage, level, T);
    }
  }

  /**
   * Generates a regular square root grid of level level with boundaries
   *
   *
   * @param storage Hashmap, that stores the grid points
   * @param level maximum level of the square root  grid (non-negative value)
   */
  void squareRoot(GridStorage& storage, level_t level) {
    if (storage.getSize() > 0) {
      throw generation_exception("storage not empty");
    }

    GridPoint point(storage.getDimension());

    for (size_t d = 0; d < storage.getDimension(); d++) {
      point.push(d, 0, 0, false);
    }

    /**
     * Change here to the following code to take the [n/2]+1 grid as small level for odd numbers(and
     * also change FullGridSet getSquare method)
     * int small_level=ceil(level/2);
     * if (level%2==0) level--;
     * */
    int small_level = level / 2;
    this->square_rec(storage, point, storage.getDimension() - 1, level, small_level, false, 0);
  }
  /**
   * Generates a truncated boundary grid containing all gridpoints with li<l-k and |l|<l+(dim-1)*k
   *
   *
   * @param storage Hashmap, that stores the grid points
   * @param level maximum level of the square root  grid (non-negative value)
   * @param k the parameter which determines the maximum level of the gridpoints for every dimension
   */
  void truncated(GridStorage& storage, level_t level, level_t k) {
    if (storage.getSize() > 0) {
      throw generation_exception("storage not empty");
    }

    GridPoint point(storage.getDimension());

    for (size_t d = 0; d < storage.getDimension(); d++) {
      point.push(d, 0, 0, false);
    }

    size_t mydim = storage.getDimension();

    point.setLeaf(true);
    trunc_rec(storage, point, (mydim - 1), static_cast<level_t>(mydim) * k,
              level + k * (static_cast<level_t>(mydim) - 1), k);
  }

 protected:
  /**
   * Generate a regular sparse grid iteratively (much faster than recursively)
   * without grid points on the boundary.
   *
   * @param storage pointer to storage object into which the grid points should be stored
   * @param n level of regular sparse grid
   * @param T modifier for subgrid selection, T = 0 implies standard sparse grid.
   *        For further information see Griebel and Knapek's paper
   *        optimized tensor-product approximation spaces
   */
  void regular_iter(GridStorage& storage, level_t n, double T = 0) {
    if (storage.getDimension() == 0) return;

    GridPoint idx_1d(storage.getDimension());

    for (size_t d = 0; d < storage.getDimension(); d++) {
      idx_1d.push(d, 1, 1, false);
    }

    // Generate 1D grid in first dimension
    for (level_t l = 1; l <= n; l++) {
      for (index_t i = 1; i < static_cast<index_t>(1 << l); i += 2) {
        if (l == n) {
          idx_1d.push(0, l, i, true);
        } else {
          idx_1d.push(0, l, i, false);
        }
        storage.insert(idx_1d);
      }
    }

    // Generate grid points in all other dimensions:
    // loop dim times over intermediate grid, take all grid points and
    // modify them in current dimension d
    for (size_t d = 1; d < storage.getDimension(); d++) {
      // current size
      size_t grid_size = storage.getSize();

      // loop over all current grid points
      for (size_t g = 0; g < grid_size; g++) {
        bool first = true;
        GridPoint idx(storage.getPoint(g));

        level_t level_sum = idx.getLevelSum() - 1;
        level_t level_max = idx.getLevelMax();

        // add remaining level-index pairs in current dimension d
        for (level_t l = 1; (static_cast<double>(l + level_sum) - (T * std::max(l, level_max)) <=
                             static_cast<double>(n + storage.getDimension() - 1) - (T * n)) &&
                            (std::max(l, level_max) <= n);
             l++) {
          for (index_t i = 1; i < static_cast<index_t>(1 << l); i += 2) {
            // first grid point is updated, all others inserted
            if (first == false) {
              // is leaf?
              if ((l + level_sum) == n + storage.getDimension() - 1) {
                idx.push(d, l, i, true);
              } else {
                idx.push(d, l, i, false);
              }
              storage.insert(idx);
            } else {
              // is leaf?
              if ((l + level_sum) == n + storage.getDimension() - 1) {
                idx.push(d, l, i, true);
              } else {
                idx.push(d, l, i, false);
              }
              storage.update(idx, g);
              first = false;
            }
          }
        }
      }
    }
  }

  void decodeCoords(DataVector& coords, std::vector<bool>& result) {
    for (size_t i = 0; i < coords.getSize(); ++i) {
      result[i] = coords[i] != 0.5;
    }
  }

  void regular_inter_iter(GridStorage& storage, level_t n,
                          const std::unordered_set<std::vector<bool>>& terms, double T = 0) {
    if (storage.getDimension() == 0) return;

    auto coords = DataVector(storage.getDimension());
    auto coordsBool = std::vector<bool>(storage.getDimension());

    GridPoint idx_1d(storage.getDimension());

    for (size_t d = 0; d < storage.getDimension(); d++) {
      idx_1d.push(d, 1, 1, false);
    }

    // Generate 1D grid in first dimension
    for (level_t l = 1; l <= n; l++) {
      for (index_t i = 1; i < static_cast<index_t>(1 << l); i += 2) {
        if (l == n) {
          idx_1d.push(0, l, i, true);
        } else {
          idx_1d.push(0, l, i, false);
        }
        storage.insert(idx_1d);
      }
    }

    // Generate grid points in all other dimensions:
    // loop dim times over intermediate grid, take all grid points and
    // modify them in current dimension d
    for (size_t d = 1; d < storage.getDimension(); d++) {
      // current size
      size_t grid_size = storage.getSize();

      // loop over all current grid points
      for (size_t g = 0; g < grid_size; g++) {
        bool first = true;
        GridPoint idx(storage.getPoint(g));

        level_t level_sum = idx.getLevelSum() - 1;
        level_t level_max = idx.getLevelMax();

        // add remaining level-index pairs in current dimension d
        for (level_t l = 1; (static_cast<double>(l + level_sum) - (T * std::max(l, level_max)) <=
                             static_cast<double>(n + storage.getDimension() - 1) - (T * n)) &&
                            (std::max(l, level_max) <= n);
             l++) {
          for (index_t i = 1; i < static_cast<index_t>(1 << l); i += 2) {
            // first grid point is updated, all others inserted
            if (first == false) {
              // is leaf?
              if ((l + level_sum) == n + storage.getDimension() - 1) {
                idx.push(d, l, i, true);
              } else {
                idx.push(d, l, i, false);
              }
              idx.rehash();
              idx.getStandardCoordinates(coords);
              decodeCoords(coords, coordsBool);
              if (terms.find(coordsBool) != terms.end()) {
                storage.insert(idx);
              }
            } else {
              // is leaf?
              if ((l + level_sum) == n + storage.getDimension() - 1) {
                idx.push(d, l, i, true);
              } else {
                idx.push(d, l, i, false);
              }
              idx.getStandardCoordinates(coords);
              decodeCoords(coords, coordsBool);
              if (terms.find(coordsBool) != terms.end()) {
                storage.update(idx, g);
                first = false;
              }
            }
          }
        }
      }
    }
  }

  void cliques_iter(GridStorage& storage, level_t n, size_t clique_size, double T = 0) {
    if (storage.getDimension() == 0) return;

    GridPoint idx_1d(storage.getDimension());

    for (size_t d = 0; d < storage.getDimension(); d++) {
      idx_1d.push(d, 1, 1, false);
    }

    // Generate 1D grid in first dimension
    for (level_t l = 1; l <= n; l++) {
      for (index_t i = 1; i < static_cast<index_t>(1 << l); i += 2) {
        if (l == n) {
          idx_1d.push(0, l, i, true);
        } else {
          idx_1d.push(0, l, i, false);
        }

        storage.insert(idx_1d);
      }
    }

    // Generate grid points in all other dimensions:
    // loop dim times over intermediate grid,
    // take all grid points and modify them in current dimension d
    for (size_t d = 1; d < storage.getDimension(); d++) {
      // current size
      size_t grid_size = storage.getSize();
      size_t clique_num = d / clique_size;

      // loop over all current grid points
      for (size_t g = 0; g < grid_size; g++) {
        bool first = true;
        bool skip = false;
        GridPoint idx(storage.getPoint(g));

        // calculate current level-sum - 1
        level_t level_sum = idx.getLevelSum() - 1;

        for (size_t dt = 0; dt < clique_size * clique_num && dt < d; dt++) {
          // if the level in dt dimension > 1, ignore the point and continue
          level_t lt;
          index_t it;  // level and index in the particular dimension
          idx.get(dt, lt, it);

          if (lt > 1) {
            skip = true;
            break;
          }
        }

        if (skip) {
          continue;
        }

        level_t level_max = idx.getLevelMax();
        // add remaining level-index pairs in current dimension d
        // as mentioned before T adjusts the granularity of the grid
        for (level_t l = 1; (static_cast<double>(l + level_sum) - (T * std::max(l, level_max)) <=
                             static_cast<double>(n + storage.getDimension() - 1) - (T * n)) &&
                            (std::max(l, level_max) <= n);
             l++) {
          for (index_t i = 1; i < static_cast<index_t>(1 << l); i += 2) {
            // first grid point is updated, all others inserted
            if (first == false) {
              // is leaf?
              if ((l + level_sum) == n + storage.getDimension() - 1) {
                idx.push(d, l, i, true);
              } else {
                idx.push(d, l, i, false);
              }

              storage.insert(idx);
            } else {
              // is leaf?
              if ((l + level_sum) == n + storage.getDimension() - 1) {
                idx.push(d, l, i, true);
              } else {
                idx.push(d, l, i, false);
              }

              storage.update(idx, g);
              first = false;
            }
          }
        }
      }
    }
  }

  /**
   * Generate a regular sparse grid iteratively (much faster than
         * recursively) with truncated boundary, i.e.,
         * the sparse grid on the \f$(d-1)\f$-dimensional faces of
         * \f$[0, 1]^d\f$ has a coarser level than the main axes
         * \f$x_t = 0.5, t = 1, ..., d\f$.
         *
         * The function adds all hierarchical subspaces \f$W_{\vec{\ell}}\f$
         * where
         * * \f$\norm{\vec{\ell}}_1 \le n + d - 1\f$ with
         *   \f$\forall_t\; \ell_t \ge 1\f$,
         * * \f$\norm{\vec{\ell}}_1 \le n + d - b - N_{\vec{\ell}}\f$ with
         *   \f$N_{\vec{\ell}} := |\{t \mid \ell_t = 0\}| \ge 1\f$, or
         * * \f$\vec{\ell} = \vec{0}\f$.
         *
         * The previous implementation inserted the 1D boundary grid points
         * at higher levels (e.g., at boundaryLevel = 2), which
         * led to the effect that corner points were missing in
         * higher-dimensional grids: For example, if \f$d = 2\f$,
         * \f$n = 3\f$, and \f$\mathtt{boundaryLevel} = 3\f$,
         * then the four corners had level sum
         * \f$2 \cdot \mathtt{boundaryLevel} = 6\f$ (which is greater than
         * \f$n + d - 1 = 4\f$), thus they were missing in the sparse grid.
         * The midpoints of the four edges, however, had level sum
         * \f$\mathtt{boundaryLevel} + 1 = 4 \le n + d - 1\f$
         * and were thus included in the grid. To get the corners into the
         * grid, too, one would have to choose \f$n \ge 4\f$.
         *
         * In contrast, the new implementation makes sure that the corners
         * will always appear first in the grid when increasing the level
         * \f$n = 1, 2, 3, \dotsc\f$ of the regular grid.
   *
   * @param storage       pointer to storage object into which
   *                      the grid points should be stored
   * @param n             level of regular sparse grid
   * @param boundaryLevel 1 + how much levels the boundary is coarser than
   *                      the main axes, 1 means same level,
   *                      2 means one level coarser, etc.; must be >= 1
   * @param T             modifier for subgrid selection, T = 0 implies standard sparse grid.
   *                      For further information see Griebel and Knapek's paper
   *                      optimized tensor-product approximation spaces
   */
  void regular_boundary_truncated_iter(GridStorage& storage, level_t n, level_t boundaryLevel = 1,
                                       double T = 0) {
    const size_t dim = storage.getDimension();

    if (dim == 0) {
      return;
    }

    GridPoint idx_1d(dim);

    for (size_t d = 0; d < dim; d++) {
      idx_1d.push(d, 1, 1, false);
    }

    // generate boundary basis functions
    {
      idx_1d.push(0, 0, 0, false);
      storage.insert(idx_1d);

      idx_1d.push(0, 0, 1, false);
      storage.insert(idx_1d);
    }

    // generate 1D grid in first dimension
    for (level_t l = 1; l <= n; l++) {
      // generate inner basis function
      for (index_t i = 1; i < static_cast<index_t>(1) << l; i += 2) {
        if (l == n) {
          idx_1d.push(0, l, i, true);
        } else {
          idx_1d.push(0, l, i, false);
        }

        storage.insert(idx_1d);
      }
    }

    // Generate grid points in all other dimensions:
    // loop dim times over intermediate grid,
    // take all grid points and modify them in current dimension d
    for (size_t d = 1; d < dim; d++) {
      // current size
      const size_t gridSize = storage.getSize();
      // curDim is new dimension of the grid points to be inserted
      const level_t curDim = static_cast<level_t>(d + 1);

      // loop over all current grid points
      for (size_t g = 0; g < gridSize; g++) {
        level_t levelSum = 0;
        level_t numberOfZeroLevels = 0;
        GridPoint idx(storage.getPoint(g));
        bool firstPoint = true;

        // calculate level sum and count number of zero levels
        for (size_t sd = 0; sd < d; sd++) {
          level_t tmp = idx.getLevel(sd);

          if (tmp == 0) {
            numberOfZeroLevels++;
          }

          levelSum += tmp;
        }

        // generate boundary basis functions,
        // but only if levelSum <=
        // n + curDim - boundaryLevel - (numberOfZeroLevels + 1)
        // (the +1 comes from the fact that the newly generated functions
        // will have an additional zero in the d-th dimension)
        if ((levelSum + boundaryLevel + numberOfZeroLevels + 1 <= n + curDim) ||
            (numberOfZeroLevels == curDim - 1)) {
          idx.push(d, 0, 0, false);
          storage.update(idx, g);

          idx.push(d, 0, 1, false);
          storage.insert(idx);

          firstPoint = false;
        }

        double upperBound;

        // choose upper bound of level sum according whether
        // the new basis function is an interior or a boundary function
        // (the loop below skips l = 0 as the boundary points
        // have been inserted a few lines above)
        if (numberOfZeroLevels > 0) {
          // check if upperBound would be negative
          // (we're working with unsigned integers here)
          if (n + curDim < boundaryLevel + numberOfZeroLevels) {
            continue;
          } else {
            // upper bound for boundary basis functions
            upperBound = static_cast<double>(n + curDim - numberOfZeroLevels - boundaryLevel);
          }
        } else {
          // upper bound for interior basis functions
          upperBound = static_cast<double>(n + curDim - 1);
        }
        upperBound -= T * n;
        level_t level_max = idx.getLevelMax();

        for (level_t l = 1;
             (static_cast<double>(l + levelSum) - (T * std::max(l, level_max)) <= upperBound) &&
             (std::max(l, level_max) <= n);
             l++) {
          // generate inner basis functions
          for (index_t i = 1; i < static_cast<index_t>(1) << l; i += 2) {
            if ((l + levelSum) == n + dim - 1) {
              idx.push(d, l, i, (numberOfZeroLevels == 0));
            } else {
              idx.push(d, l, i, false);
            }

            if (firstPoint) {
              storage.update(idx, g);
              firstPoint = false;
            } else {
              storage.insert(idx);
            }
          }
        }
      }
    }
  }

  /**
   * Generate a regular sparse grid iteratively (much faster than recursively)
   * with periodic boundary.
   *
   * @param storage Pointer to storage object into which the grid points should be stored
   * @param n Level of regular sparse grid
   * @param T modifier for subgrid selection, T = 0 implies standard sparse grid.
   *           For further information see Griebel and Knapek's paper
   *           optimized tensor-product approximation spaces
   */
  void regular_periodic_boundary_iter(GridStorage& storage, level_t n, double T = 0) {
    if (storage.getDimension() == 0) return;

    GridPoint idx_1d(storage.getDimension());

    for (size_t d = 0; d < storage.getDimension(); d++) {
      idx_1d.push(d, 1, 1, false);
    }

    // Generate 1D grid in first dimension
    for (level_t l = 1; l <= n; l++) {
      // generate boundary basis functions
      if (l == 1) {
        idx_1d.push(0, 0, 0, false);
        storage.insert(idx_1d);
      }

      for (index_t i = 1; i < static_cast<index_t>(1 << l); i += 2) {
        if (l == n) {
          idx_1d.push(0, l, i, true);
        } else {
          idx_1d.push(0, l, i, false);
        }

        storage.insert(idx_1d);
      }
    }

    // Generate grid points in all other dimensions:
    // loop dim times over intermediate grid, take all grid points and
    // modify them in current dimension d
    for (size_t d = 1; d < storage.getDimension(); d++) {
      // current size
      size_t grid_size = storage.getSize();

      // loop over all current grid points
      for (size_t g = 0; g < grid_size; g++) {
        bool first = true;
        GridPoint idx(storage.getPoint(g));

        // Calculate level-sum
        level_t level_sum = idx.getLevelSum() - 1;

        for (size_t sd = 0; sd < d; sd++) {
          if (idx.getLevel(sd) == 0) level_sum += 1;
        }

        level_t level_max = idx.getLevelMax();
        // add remaining level-index pairs in current dimension d
        // as mentioned before T adjusts the granularity of the grid
        for (level_t l = 1; (static_cast<double>(l + level_sum) - (T * std::max(l, level_max)) <=
                             static_cast<double>(n + storage.getDimension() - 1) - (T * n)) &&
                            (std::max(l, level_max) <= n);
             l++) {
          if (l == 1) {
            idx.push(d, 0, 0, false);
            storage.insert(idx);
          }

          for (index_t i = 1; i < static_cast<index_t>(1 << l); i += 2) {
            // first grid point is updated, all others inserted
            if (first == false) {
              // is leaf?
              if ((l + level_sum) == n + storage.getDimension() - 1) {
                idx.push(d, l, i, true);
              } else {
                idx.push(d, l, i, false);
              }

              storage.insert(idx);
            } else {
              // is leaf?
              if ((l + level_sum) == n + storage.getDimension() - 1) {
                idx.push(d, l, i, true);
              } else {
                idx.push(d, l, i, false);
              }

              storage.update(idx, g);
              first = false;
            }
          }
        }
      }
    }
  }

  /**
   * Generate a full grid iteratively (much faster than recursively)
   * without grid points on the boundary.
   *
   * @param storage Pointer to the storage object into which the grid points should be stored
   * @param n Level of full grid
   */
  void createFullGridIterative(GridStorage& storage, level_t n) {
    if (storage.getDimension() == 0) return;

    GridPoint idx_1d(storage.getDimension());

    for (size_t d = 0; d < storage.getDimension(); d++) {
      idx_1d.push(d, 1, 1, false);
    }

    // Generate 1D grid in first dimension
    for (level_t l = 1; l <= n; l++) {
      for (index_t i = 1; i < static_cast<index_t>(1 << l); i += 2) {
        if (l == n) {
          idx_1d.push(0, l, i, true);
        } else {
          idx_1d.push(0, l, i, false);
        }

        storage.insert(idx_1d);
      }
    }

    // Generate grid points in all other dimensions:
    // loop dim times over intermediate grid, take all grid points and
    // modify them in current dimension d
    for (size_t d = 1; d < storage.getDimension(); d++) {
      // current size
      size_t grid_size = storage.getSize();

      // loop over all current grid points
      for (size_t g = 0; g < grid_size; g++) {
        bool first = true;
        GridPoint idx(storage.getPoint(g));

        // add remaining level-index pairs in current dimension d
        for (level_t l = 1; l <= n; l++) {
          // for leaf check, set current level to l
          idx.push(d, l, 1, false);

          for (index_t i = 1; i < static_cast<index_t>(1 << l); i += 2) {
            // first grid point is updated, all others inserted
            if (first == false) {
              // is leaf?
              if (idx.getLevelSum() == n * storage.getDimension()) {
                idx.push(d, l, i, true);
              } else {
                idx.push(d, l, i, false);
              }

              storage.insert(idx);
            } else {
              // is leaf?
              if (idx.getLevelSum() == n * storage.getDimension()) {
                idx.push(d, l, i, true);
              } else {
                idx.push(d, l, i, false);
              }

              storage.update(idx, g);
              first = false;
            }
          }
        }
      }
    }
  }

  void createAnisotropicFullGrid(GridStorage& storage, std::vector<size_t> v) {
    if (storage.getDimension() != v.size()) {
      throw sgpp::base::generation_exception(
          "HashGenerator:: createAnisotropicFullGrid() Storage size doesn't fit vector size");
      return;
    }

    if (storage.getDimension() == 0) {
      std::cout << "Storage dimension is zero!";
      return;
    }

    GridPoint idx_1d(storage.getDimension());

    for (size_t d = 0; d < storage.getDimension(); d++) {
      idx_1d.push(d, 1, 1, false);
    }

    // Generate 1D grid in first dimension
    for (level_t l = 1; l <= v.at(0); l++) {
      for (index_t i = 1; i < static_cast<index_t>(1 << l); i += 2) {
        if (l == v.at(0)) {
          idx_1d.push(0, l, i, true);
        } else {
          idx_1d.push(0, l, i, false);
        }

        storage.insert(idx_1d);
      }
    }

    // Generate grid points in all other dimensions:
    // loop dim times over intermediate grid, take all grid points and
    // modify them in current dimension d
    for (size_t d = 1; d < storage.getDimension(); d++) {
      // current size
      size_t grid_size = storage.getSize();

      // loop over all current grid points
      for (size_t g = 0; g < grid_size; g++) {
        bool first = true;
        GridPoint idx(storage.getPoint(g));

        // add remaining level-index pairs in current dimension d
        for (level_t l = 1; l <= v.at(d); l++) {
          // for leaf check, set current level to l
          idx.push(d, l, 1, false);

          for (index_t i = 1; i < static_cast<index_t>(1 << l); i += 2) {
            // first grid point is updated, all others inserted
            if (first == false) {
              // is leaf?
              if (idx.getLevelSum() == v.at(d) * storage.getDimension()) {
                idx.push(d, l, i, true);
              } else {
                idx.push(d, l, i, false);
              }

              storage.insert(idx);
            } else {
              // is leaf?
              if (idx.getLevelSum() == v.at(d) * storage.getDimension()) {
                idx.push(d, l, i, true);
              } else {
                idx.push(d, l, i, false);
              }

              storage.update(idx, g);
              first = false;
            }
          }
        }
      }
    }
  }

  /**
   * Generate a full grid iteratively (much faster than recursively)
   * with truncated boundary.
   *
   * @param storage Pointer to the storage object into which the grid points should be stored
   * @param n Level of full grid
   */
  void createFullGridTruncatedIterative(GridStorage& storage, level_t n) {
    if (storage.getDimension() == 0) return;

    GridPoint idx_1d(storage.getDimension());

    for (size_t d = 0; d < storage.getDimension(); d++) {
      idx_1d.push(d, 1, 1, false);
    }

    // Generate 1D grid in first dimension
    for (level_t l = 1; l <= n; l++) {
      // generate boundary basis functions
      if (l == 1) {
        idx_1d.push(0, 0, 0, false);
        storage.insert(idx_1d);
        idx_1d.push(0, 0, 1, false);
        storage.insert(idx_1d);
      }

      // generate inner basis function
      for (index_t i = 1; i < static_cast<index_t>(1 << l); i += 2) {
        if (l == n) {
          idx_1d.push(0, l, i, true);
        } else {
          idx_1d.push(0, l, i, false);
        }

        storage.insert(idx_1d);
      }
    }

    // Generate grid points in all other dimensions:
    // loop dim times over intermediate grid, take all grid points and
    // modify them in current dimension d
    for (size_t d = 1; d < storage.getDimension(); d++) {
      // current size
      size_t grid_size = storage.getSize();

      // loop over all current grid points
      for (size_t g = 0; g < grid_size; g++) {
        GridPoint idx(storage.getPoint(g));

        // add remaining level-index pairs in current dimension d
        for (level_t l = 1; l <= n; l++) {
          // generate boundary basis functions
          if (l == 1) {
            idx.push(d, 0, 0, false);
            storage.update(idx, g);

            idx.push(d, 0, 1, false);
            storage.insert(idx);
          }

          // for leaf check, set current level to l
          idx.push(d, l, 1, false);

          // generate inner basis functions
          for (index_t i = 1; i < static_cast<index_t>(1 << l); i += 2) {
            // is leaf?
            if (idx.getLevelSum() == n * storage.getDimension()) {
              idx.push(d, l, i, true);
            } else {
              idx.push(d, l, i, false);
            }

            storage.insert(idx);
          }
        }
      }
    }
  }

  //  /**
  //   * recursive construction of the spare grid without boundaries
  //   *
  //   * @param storage hashmap that stores the grid points
  //   * @param index point's index
  //   * @param current_dim current working dimension
  //   * @param current_level current level in this construction step
  //   * @param level maximum level of the sparse grid
  //   */
  //  void regular_rec(GridStorage& storage, GridPoint& index,
  //  size_t current_dim, level_t current_level, level_t level)
  //  {
  //    if(current_dim == 0)
  //    {
  //      regular_rec_1d(storage, index, current_level, level);
  //    }
  //    else
  //    {
  //      index_t source_index;
  //      level_t source_level;
  //
  //      index.get(current_dim, source_level, source_index);
  //
  //      if(current_level <= level)
  //      {
  //        // set Leaf option of index
  //        if (current_level == level)
  //        {
  //          index.setLeaf(true);
  //        }
  //        else
  //        {
  //          index.setLeaf(false);
  //        }
  //
  //        // d-1 recursion
  //        this->regular_rec(storage, index, current_dim - 1, current_level,
  //        level);
  //      }
  //
  //      if(current_level < level)
  //      {
  //        // current_level + 1 recursion
  //        index.push(current_dim, source_level + 1, 2*source_index - 1);
  //        this->regular_rec(storage, index, current_dim, current_level + 1,
  //        level);
  //
  //        index.push(current_dim, source_level + 1, 2*source_index + 1);
  //        this->regular_rec(storage, index, current_dim, current_level + 1,
  //        level);
  //      }
  //
  //      index.push(current_dim, source_level, source_index);
  //    }
  //  }
  //
  //  /**
  //   * generate points of the last dimension (dim == 0), without boundaries
  //   *
  //   * @param storage the hashmap that stores the grid points
  //   * @param index point's index that should be created on the grid
  //   * @param current_level current level of the grid generation
  //   * @param level maximum level of grid
  //   */
  //  void regular_rec_1d(GridStorage& storage, GridPoint& index,
  //  level_t current_level, level_t level)
  //  {
  //        for(level_t l = 1; l <= level-current_level + 1; l++)
  //        {
  //            if (l == level-current_level+1)
  //            {
  //              for(index_t i = 1;
  //              static_cast<int>(i) <= static_cast<int>(1<<(l-1)); i++)
  //                {
  //                    index.push(0, l, 2*i-1, true);
  //                    storage.insert(index);
  //                }
  //            }
  //            else
  //            {
  //              for(index_t i = 1;
  //              static_cast<int>(i) <= static_cast<int>(1<<(l-1)); i++)
  //                {
  //                    index.push(0, l, 2*i-1, false);
  //                    storage.insert(index);
  //                }
  //            }
  //        }
  //  }

  /**
   * recursive construction of the spare grid with boundaries, pentagon cut
   *
   * @param storage hashmap that stores the grid points
   * @param index point's index
   * @param current_dim current working dimension
   * @param current_level current level in this construction step
   * @param level maximum level of the sparse grid
   * @param bLevelZero specifies if the current index has a level zero component
   */

  void boundaries_truncated_rec(GridStorage& storage, GridPoint& index, size_t current_dim,
                                level_t current_level, level_t level, bool bLevelZero) {
    if (current_dim == 0) {
      boundaries_Truncated_rec_1d(storage, index, current_level, level, bLevelZero);
    } else {
      index_t source_index;
      level_t source_level;

      index.get(current_dim, source_level, source_index);

      if (current_level <= level) {
        // set Leaf option of index
        bool bLeafProperty = false;

        if (current_level == level) {
          bLeafProperty = true;
        } else {
          bLeafProperty = false;
        }

        if (source_level == 1) {
          index.push(current_dim, 0, 0, false);
          this->boundaries_truncated_rec(storage, index, current_dim - 1, current_level, level,
                                         true);

          index.push(current_dim, 0, 1, false);
          this->boundaries_truncated_rec(storage, index, current_dim - 1, current_level, level,
                                         true);

          index.push(current_dim, source_level, source_index);
        }

        // d-1 recursion
        index.setLeaf(bLeafProperty);
        this->boundaries_truncated_rec(storage, index, current_dim - 1, current_level, level,
                                       bLevelZero);
      }

      if (current_level < level) {
        index.push(current_dim, source_level + 1, 2 * source_index - 1);
        this->boundaries_truncated_rec(storage, index, current_dim, current_level + 1, level,
                                       bLevelZero);

        index.push(current_dim, source_level + 1, 2 * source_index + 1);
        this->boundaries_truncated_rec(storage, index, current_dim, current_level + 1, level,
                                       bLevelZero);
      }

      index.push(current_dim, source_level, source_index);
    }
  }

  /**
   * generate points of the last dimension (dim == 0), version of pentagon cut in
   * sub space scheme
   *
   * @param storage the hashmap that stores the grid points
   * @param index point's index that should be created on the grid
   * @param current_level current level of the grid generation
   * @param level maximum level of grid
   * @param bLevelZero specifies if the current index has a level zero component
   */

  void boundaries_Truncated_rec_1d(GridStorage& storage, GridPoint& index, level_t current_level,
                                   level_t level, bool bLevelZero) {
    bool bLevelGreaterZero = !bLevelZero;

    for (level_t l = 0; l <= level - current_level + 1; l++) {
      if (l == level - current_level + 1) {
        if (l == 0) {
          index.push(0, 0, 0, false);
          storage.insert(index);
          index.push(0, 0, 1, false);
          storage.insert(index);
        } else {
          for (index_t i = 1; static_cast<int>(i) <= static_cast<int>(1 << (l - 1)); i++) {
            index.push(0, l, 2 * i - 1, (true && bLevelGreaterZero));
            storage.insert(index);
          }
        }
      } else {
        if (l == 0) {
          index.push(0, 0, 0, false);
          storage.insert(index);
          index.push(0, 0, 1, false);
          storage.insert(index);
        } else {
          for (index_t i = 1; static_cast<int>(i) <= static_cast<int>(1 << (l - 1)); i++) {
            index.push(0, l, 2 * i - 1, false);
            storage.insert(index);
          }
        }
      }
    }
  }

  /**
   * recursive construction of the spare grid with boundaries, classic level 0 approach, only for
   *level 0 and 1
   *
   * @param storage hashmap that stores the grid points
   * @param index point's index
   * @param current_dim current working dimension
   * @param current_level current level in this construction step
   * @param level maximum level of the sparse grid
   */

  void boundaries_rec(GridStorage& storage, GridPoint& index, size_t current_dim,
                      level_t current_level, level_t level) {
    index_t source_index;
    level_t source_level;

    index.get(current_dim, source_level, source_index);

    if (current_level <= level) {
      // set Leaf option of index
      bool bSaveLeafProperty = index.isLeaf();
      bool bLeafProperty = false;

      if (current_level == level) {
        bLeafProperty = true;
      } else {
        bLeafProperty = false;
      }

      // d-1 recursion
      if (source_level == 0) {
        if (current_dim == 0) {
          index.push(0, 0, 0, bLeafProperty);
          storage.insert(index);
          index.push(0, 0, 1, bLeafProperty);
          storage.insert(index);

          index.push(current_dim, source_level, source_index, bSaveLeafProperty);
        } else {
          index.push(current_dim, 0, 0, bLeafProperty);
          this->boundaries_rec(storage, index, current_dim - 1, current_level, level);

          index.push(current_dim, 0, 1, bLeafProperty);
          this->boundaries_rec(storage, index, current_dim - 1, current_level, level);

          index.push(current_dim, source_level, source_index, bSaveLeafProperty);
        }
      } else {
        index.setLeaf(bLeafProperty);

        if (current_dim == 0) {
          storage.insert(index);
        } else {
          this->boundaries_rec(storage, index, current_dim - 1, current_level, level);
        }

        index.setLeaf(bSaveLeafProperty);
      }
    }

    if (current_level < level) {
      if (source_level == 0 && source_index == 0) {
        index.push(current_dim, source_level + 1, 1);
        this->boundaries_rec(storage, index, current_dim, current_level + 1, level);
      } else {
        index.push(current_dim, source_level + 1, 2 * source_index - 1);
        this->boundaries_rec(storage, index, current_dim, current_level + 1, level);

        index.push(current_dim, source_level + 1, 2 * source_index + 1);
        this->boundaries_rec(storage, index, current_dim, current_level + 1, level);
      }
    }

    index.push(current_dim, source_level, source_index);
  }
  /**
   * recursive construction of a square root grid with boundaries
   *
   * @param storage hashmap that stores the grid points
   * @param index point's index
   * @param current_dim current working dimension
   * @param level maximum level of the square root grid
   * @param small_level level of coarsest descretization
   * @param tail true if there is a level of the index>level/2
   * @param sum sum of all levels
   */

  void square_rec(GridStorage& storage, GridPoint& index, size_t current_dim, level_t level,
                  level_t small_level, bool tail, size_t sum) {
    index_t source_index;
    level_t source_level;

    index.get(current_dim, source_level, source_index);
    bool newtail = tail || (source_level > small_level);

    // d-1 recursion
    if (source_level == 0) {
      if (current_dim == 0) {
        index.setLeaf(false);
        index.push(0, 0, 0);
        storage.insert(index);
        index.push(0, 0, 1);
        storage.insert(index);

        index.push(current_dim, source_level, source_index);
      } else {
        // we want all nodes,
        // so we initialize the unassigned levels to start from 0
        for (size_t d = 0; d < current_dim; d++) index.push(d, 0, 0);

        index.push(current_dim, 0, 0);
        this->square_rec(storage, index, current_dim - 1, level, small_level, newtail, sum);

        for (size_t d = 0; d < current_dim; d++) index.push(d, 0, 0);

        index.push(current_dim, 0, 1);

        this->square_rec(storage, index, current_dim - 1, level, small_level, newtail, sum);

        index.push(current_dim, source_level, source_index);
      }
    } else {
      if (current_dim == 0) {
        /**
         * If a level of the node equals level, and all the others equal small_level, the node is a
         * leaf
         * This is equivalent to saying the sum of levels equals small_level*(dim-1)+level
         * */
        if (sum == (small_level * (storage.getDimension() - 1) + level)) {
          index.setLeaf(true);
        } else {
          index.setLeaf(false);
        }

        storage.insert(index);
      } else {
        for (size_t d = 0; d < current_dim; d++) index.push(d, 0, 0);

        this->square_rec(storage, index, current_dim - 1, level, small_level, newtail, sum);
      }
    }

    /**If the level of the node is smaller than small_level or we didn't have yet a level greater
    than small_level(!tail) and the level is smaller then level then
    we can then proceed to the next level on this dimension, otherwise we reached the maximum
    possible level*
    */
    if ((source_level < small_level) || ((!tail) && (source_level < level))) {
      if (source_level == 0 && source_index == 0) {
        index.push(current_dim, source_level + 1, 1);
        this->square_rec(storage, index, current_dim, level, small_level, tail, sum + 1);
      } else {
        index.push(current_dim, source_level + 1, 2 * source_index - 1);
        this->square_rec(storage, index, current_dim, level, small_level, tail, sum + 1);

        index.push(current_dim, source_level + 1, 2 * source_index + 1);
        this->square_rec(storage, index, current_dim, level, small_level, tail, sum + 1);
      }
    }
  }
  /**
   * recursive construction of a super truncated grid with boundaries
   *
   * @param storage hashmap that stores the grid points
   * @param index point's index
   * @param current_dim current working dimension
   * @param current_level the current level of the gridpoint so far, starts from minlevel*dim
   * @param level the maximum level of the gridpoint
   * @param minlevel the level limit given by the user(tells us which fullgrids won't be present in
   *the construction of the sparse grid)
   */

  void trunc_rec(GridStorage& storage, GridPoint& index, size_t current_dim, level_t current_level,
                 level_t level, level_t minlevel) {
    index_t source_index;
    level_t source_level;

    index.get(current_dim, source_level, source_index);
    bool bSaveLeafProperty = index.isLeaf();

    if (current_level <= level) {
      // set Leaf option of index

      bool bLeafProperty = bSaveLeafProperty;

      if (source_level < minlevel) {
        bLeafProperty = false;
      }

      // d-1 recursion
      if (source_level == 0) {
        if (current_dim == 0) {
          if (current_level < level) bLeafProperty = false;

          index.push(0, 0, 0, bLeafProperty);
          storage.insert(index);
          index.push(0, 0, 1, bLeafProperty);
          storage.insert(index);

          index.push(current_dim, source_level, source_index, bSaveLeafProperty);
        } else {
          index.push(current_dim, 0, 0, bLeafProperty);
          trunc_rec(storage, index, current_dim - 1, current_level, level, minlevel);

          index.push(current_dim, 0, 1, bLeafProperty);
          trunc_rec(storage, index, current_dim - 1, current_level, level, minlevel);

          index.push(current_dim, source_level, source_index, bSaveLeafProperty);
        }
      } else {
        if (current_dim == 0) {
          if (current_level < level) bLeafProperty = false;

          index.setLeaf(bLeafProperty);
          storage.insert(index);
        } else {
          index.setLeaf(bLeafProperty);
          trunc_rec(storage, index, current_dim - 1, current_level, level, minlevel);
        }

        index.setLeaf(bSaveLeafProperty);
      }
    }

    if (source_level < minlevel) {
      /**
       * If the source level of the node is smaller than minlevel we don't increase the variable
       * current_level(since we started with minlevel*dim)
       * This trick makes it possible to introduce all nodes with source_level&lt;minlevel without a
       * separate treatment
       * */
      if (source_level == 0 && source_index == 0) {
        index.push(current_dim, source_level + 1, 1);
        trunc_rec(storage, index, current_dim, current_level, level, minlevel);
      } else {
        index.push(current_dim, source_level + 1, 2 * source_index - 1);
        trunc_rec(storage, index, current_dim, current_level, level, minlevel);

        index.push(current_dim, source_level + 1, 2 * source_index + 1);
        trunc_rec(storage, index, current_dim, current_level, level, minlevel);
      }

    } else if (current_level < level) {
      /**
       * if the source_level is already >=minlevel we can proceed naturally and increase the
       * current_level which represents the sum of levels so far
       * */
      if (source_level == 0 && source_index == 0) {
        index.push(current_dim, source_level + 1, 1);
        trunc_rec(storage, index, current_dim, current_level + 1, level, minlevel);
      } else {
        index.push(current_dim, source_level + 1, 2 * source_index - 1);
        trunc_rec(storage, index, current_dim, current_level + 1, level, minlevel);

        index.push(current_dim, source_level + 1, 2 * source_index + 1);
        trunc_rec(storage, index, current_dim, current_level + 1, level, minlevel);
      }
    }

    index.push(current_dim, source_level, source_index, bSaveLeafProperty);
  }
};

}  // namespace base
}  // namespace sgpp

#endif /* HASHGENERATOR_HPP */
