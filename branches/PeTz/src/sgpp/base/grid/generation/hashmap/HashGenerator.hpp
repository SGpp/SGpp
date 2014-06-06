/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Dirk Pflueger (pflueged@in.tum.de), Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef HASHGENERATOR_HPP
#define HASHGENERATOR_HPP

#include "base/grid/GridStorage.hpp"

#include "base/exception/generation_exception.hpp"

#include <vector>
#include <cmath>
#include <iostream>

namespace sg {
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
     * @todo (heinecke, nice) add picture here
     *
     * 2. a modified boundary grid with level 0 and a pentagon cut
     * trough the sub space scheme.
     *
     * @todo (heinecke, nice) add picture here
     *
     * Furthermore, the creation of full grids (in the hierarchical basis) is supported.
     */
    class HashGenerator {
      public:

        typedef GridStorage::index_type index_type;
        typedef index_type::index_type index_t;
        typedef index_type::level_type level_t;

        /**
         * Generates a regular sparse grid of level levels, without boundaries
         *
         * @todo (blank) level should be of type level_t but swig doesnt want that
         *
         * @param storage Hashmap that stores the grid points
         * @param level Grid level (non-negative value)
         */
        void regular(GridStorage* storage, level_t level) {
          if (storage->size() > 0) {
            throw generation_exception("storage not empty");
          }

          this->regular_iter(storage, level);
        }

        /**
         * Generates a full grid of level @p level, without boundaries.
         *
         * @todo (blank) level should be of type level_t but swig doesnt want that
         *
         * @param storage Hashmap that stores the grid points
         * @param level Grid level (non-negative value)
         */
        void full(GridStorage* storage, level_t level) {
          if (storage->size() > 0) {
            throw generation_exception("storage not empty");
          }

          this->createFullGridIterative(storage, level);
        }

        /**
         * Generates a full grid of level @p level, with boundary grid points.
         *
         * @todo (blank) level should be of type level_t but swig doesnt want that
         *
         * @param storage Hashmap that stores the grid points
         * @param level Grid level (non-negative value)
         */
        void fullWithBoundary(GridStorage* storage, level_t level) {
          if (storage->size() > 0) {
            throw generation_exception("storage not empty");
          }

          this->createFullGridTrapezoidIterative(storage, level);
        }

        /**
         * Generates a regular sparse grid of level levels with boundaries
         *
         * @todo (blank) level should be of type level_t but swig doesnt want that
         *
         * @param storage Hashmap, that stores the grid points
         * @param level maximum level of the sparse grid (non-negative value)
         * @param bTrapezoidBoundaries true -> generate sparse grid with less points on the boundary, pentagon cut through subspace scheme
         */
        void regularWithBoundaries(GridStorage* storage, level_t level, bool bTrapezoidBoundaries) {
          if (storage->size() > 0) {
            throw generation_exception("storage not empty");
          }

          index_type index(storage->dim());

          if (bTrapezoidBoundaries == true) {
            if (level == 0) {
              for (size_t d = 0; d < storage->dim(); d++) {
                index.push(d, 0, 0, false);
              }

              this->boundaries_rec(storage, index, storage->dim() - 1, 0, 0);
            } else {
              //        for(size_t d = 0; d < storage->dim(); d++)
              //        {
              //          index.push(d, 1, 1, false);
              //        }
              //
              //        this->boundaries_trapezoid_rec(storage, index, storage->dim()-1, storage->dim(), level + storage->dim() - 1, false);

              this->regular_boundary_trapezoid_iter(storage, level);
            }
          } else {
            /* new grid generation
             *
             * for all level the same calculation of the level sum is implemented:
             * |l| <= n
             */
            for (size_t d = 0; d < storage->dim(); d++) {
              index.push(d, 0, 0, false);
            }

            this->boundaries_rec(storage, index, storage->dim() - 1, 0, level);
          }
        }
        /**
         * Generates a regular square root grid of level level with boundaries
         *
         *
         * @param storage Hashmap, that stores the grid points
         * @param level maximum level of the square root  grid (non-negative value)
         */
        void squareRoot(GridStorage* storage, level_t level) {
          if (storage->size() > 0) {
            throw generation_exception("storage not empty");
          }

          index_type index(storage->dim());

          for (size_t d = 0; d < storage->dim(); d++) {
            index.push(d, 0, 0, false);
          }

          /**
           * Change here to the following code to take the [n/2]+1 grid as small level for odd numbers(and also change FullGridSet getSquare method)
           * int small_level=ceil(level/2);
           * if (level%2==0) level--;
           * */
          int small_level = level / 2;
          this->square_rec(storage, index, storage->dim() - 1,  level, small_level, false, 0);

        }
        /**
         * Generates a truncated trapezoid boundary grid containing all gridpoints with li<l-k and |l|<l+(dim-1)*k
         *
         *
         * @param storage Hashmap, that stores the grid points
         * @param level maximum level of the square root  grid (non-negative value)
         * @param k the parameter which determines the maximum level of the gridpoints for every dimension
         */
        void truncated(GridStorage* storage, level_t level, level_t k) {
          if (storage->size() > 0) {
            throw generation_exception("storage not empty");
          }

          index_type index(storage->dim());

          for (size_t d = 0; d < storage->dim(); d++) {
            index.push(d, 0, 0, false);
          }

          size_t mydim = storage->dim();
          //level_t dimt = *reinterpret_cast<level_t*>(&mydim);
          index.setLeaf(true);
          trunc_rec(storage, index, (mydim - 1), *reinterpret_cast<level_t*>(&mydim) * k,
                    level + k * (*reinterpret_cast<level_t*>(&mydim) - 1), k);
        }

      protected:

        /**
         * Generate a regular sparse grid iteratively (much faster than recursively)
         * without grid points on the boundary.
         *
         * @param storage pointer to storage object into which the grid points should be stored
         * @param n level of regular sparse grid
         */
        void regular_iter(GridStorage* storage, level_t n) {
          if (storage->dim() == 0)
            return;

          index_type idx_1d(storage->dim());

          for (size_t d = 0; d < storage->dim(); d++) {
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

              storage->insert(idx_1d);
            }
          }

          // Generate grid points in all other dimensions:
          // loop dim times over intermediate grid, take all grid points and modify them
          // in current dimension d
          for (size_t d = 1; d < storage->dim(); d++) {
            // current size
            size_t grid_size = storage->size();

            // loop over all current grid points
            for (size_t g = 0; g < grid_size; g++) {
              bool first = true;
              index_type idx(storage->get(g));

              // calculate current level-sum - 1
              level_t level_sum = idx.getLevelSum() - 1;

              // add remaining level-index pairs in current dimension d
              for (level_t l = 1; l + level_sum <= n + storage->dim() - 1; l++) {
                for (index_t i = 1; i < static_cast<index_t>(1 << l); i += 2) {
                  // first grid point is updated, all others inserted
                  if (first == false) {
                    // is leaf?
                    if ((l + level_sum) == n + storage->dim() - 1) {
                      idx.push(d, l, i, true);
                    } else {
                      idx.push(d, l, i, false);
                    }

                    storage->insert(idx);
                  } else {
                    // is leaf?
                    if ((l + level_sum) == n + storage->dim() - 1) {
                      idx.push(d, l, i, true);
                    } else {
                      idx.push(d, l, i, false);
                    }

                    storage->update(idx, g);
                    first = false;
                  }
                }
              }
            }
          }
        }

        /**
         * Generate a regular sparse grid iteratively (much faster than recursively)
         * with trapezoidal boundary.
         *
         * @param storage Pointer to storage object into which the grid points should be stored
         * @param n Level of regular sparse grid
         */
        void regular_boundary_trapezoid_iter(GridStorage* storage, level_t n) {
          if (storage->dim() == 0)
            return;

          index_type idx_1d(storage->dim());

          for (size_t d = 0; d < storage->dim(); d++) {
            idx_1d.push(d, 1, 1, false);
          }

          // Generate 1D grid in first dimension
          for (level_t l = 1; l <= n; l++) {
            // generate boundary basis functions
            if (l == 1) {
              idx_1d.push(0, 0, 0, false);
              storage->insert(idx_1d);
              idx_1d.push(0, 0, 1, false);
              storage->insert(idx_1d);
            }

            // generate inner basis function
            for (index_t i = 1; i < static_cast<index_t>(1 << l); i += 2) {
              if (l == n) {
                idx_1d.push(0, l, i, true);
              } else {
                idx_1d.push(0, l, i, false);
              }

              storage->insert(idx_1d);
            }
          }

          // Generate grid points in all other dimensions:
          // loop dim times over intermediate grid, take all grid points and modify them
          // in current dimension d
          for (size_t d = 1; d < storage->dim(); d++) {
            // current size
            size_t grid_size = storage->size();

            // loop over all current grid points
            for (size_t g = 0; g < grid_size; g++) {
              // add missing Level 1
              level_t level_sum = (level_t)((storage->dim() - 1) - d);
              bool has_level_zero = false;
              index_type idx(storage->get(g));

              // Calculate level-sum
              for (size_t sd = 0; sd < d; sd++) {
                level_t tmp = idx.getLevel(sd);

                if (tmp == 0) {
                  tmp = 1;
                  has_level_zero = true;
                }

                level_sum += tmp;
              }

              for (level_t l = 1; l + level_sum <= n + storage->dim() - 1; l++) {
                // generate boundary basis functions
                if (l == 1) {
                  idx.push(d, 0, 0, false);
                  storage->update(idx, g);

                  idx.push(d, 0, 1, false);
                  storage->insert(idx);
                }

                // generate inner basis functions
                for (index_t i = 1; i < static_cast<index_t>(1 << l); i += 2) {
                  if ((l + level_sum) == n + storage->dim() - 1) {
                    if (has_level_zero == false) {
                      idx.push(d, l, i, true);
                    } else {
                      idx.push(d, l, i, false);
                    }
                  } else {
                    idx.push(d, l, i, false);
                  }

                  storage->insert(idx);
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
        void createFullGridIterative(GridStorage* storage, level_t n) {
          if (storage->dim() == 0)
            return;

          index_type idx_1d(storage->dim());

          for (size_t d = 0; d < storage->dim(); d++) {
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

              storage->insert(idx_1d);
            }
          }

          // Generate grid points in all other dimensions:
          // loop dim times over intermediate grid, take all grid points and modify them
          // in current dimension d
          for (size_t d = 1; d < storage->dim(); d++) {
            // current size
            size_t grid_size = storage->size();

            // loop over all current grid points
            for (size_t g = 0; g < grid_size; g++) {
              bool first = true;
              index_type idx(storage->get(g));

              // add remaining level-index pairs in current dimension d
              for (level_t l = 1; l <= n; l++) {
                // for leaf check, set current level to l
                idx.push(d, l, 1, false);

                for (index_t i = 1; i < static_cast<index_t>(1 << l); i += 2) {
                  // first grid point is updated, all others inserted
                  if (first == false) {
                    // is leaf?
                    if (idx.getLevelSum() == n * storage->dim()) {
                      idx.push(d, l, i, true);
                    } else {
                      idx.push(d, l, i, false);
                    }

                    storage->insert(idx);
                  } else {
                    // is leaf?
                    if (idx.getLevelSum() == n * storage->dim()) {
                      idx.push(d, l, i, true);
                    } else {
                      idx.push(d, l, i, false);
                    }

                    storage->update(idx, g);
                    first = false;
                  }
                }
              }
            }
          }
        }

        /**
         * Generate a full grid iteratively (much faster than recursively)
         * with trapezoidal boundary.
         *
         * @param storage Pointer to the storage object into which the grid points should be stored
         * @param n Level of full grid
         */
        void createFullGridTrapezoidIterative(GridStorage* storage, level_t n) {
          if (storage->dim() == 0)
            return;

          index_type idx_1d(storage->dim());

          for (size_t d = 0; d < storage->dim(); d++) {
            idx_1d.push(d, 1, 1, false);
          }

          // Generate 1D grid in first dimension
          for (level_t l = 1; l <= n; l++) {
            // generate boundary basis functions
            if (l == 1) {
              idx_1d.push(0, 0, 0, false);
              storage->insert(idx_1d);
              idx_1d.push(0, 0, 1, false);
              storage->insert(idx_1d);
            }

            // generate inner basis function
            for (index_t i = 1; i < static_cast<index_t>(1 << l); i += 2) {
              if (l == n) {
                idx_1d.push(0, l, i, true);
              } else {
                idx_1d.push(0, l, i, false);
              }

              storage->insert(idx_1d);
            }
          }

          // Generate grid points in all other dimensions:
          // loop dim times over intermediate grid, take all grid points and modify them
          // in current dimension d
          for (size_t d = 1; d < storage->dim(); d++) {
            // current size
            size_t grid_size = storage->size();

            // loop over all current grid points
            for (size_t g = 0; g < grid_size; g++) {
              index_type idx(storage->get(g));

              // add remaining level-index pairs in current dimension d
              for (level_t l = 1; l <= n; l++) {
                // generate boundary basis functions
                if (l == 1) {
                  idx.push(d, 0, 0, false);
                  storage->update(idx, g);

                  idx.push(d, 0, 1, false);
                  storage->insert(idx);
                }

                // for leaf check, set current level to l
                idx.push(d, l, 1, false);

                // generate inner basis functions
                for (index_t i = 1; i < static_cast<index_t>(1 << l); i += 2) {
                  // is leaf?
                  if (idx.getLevelSum() == n * storage->dim()) {
                    idx.push(d, l, i, true);
                  } else {
                    idx.push(d, l, i, false);
                  }

                  storage->insert(idx);
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
        //  void regular_rec(GridStorage* storage, index_type& index, size_t current_dim, level_t current_level, level_t level)
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
        //        this->regular_rec(storage, index, current_dim - 1, current_level, level);
        //      }
        //
        //      if(current_level < level)
        //      {
        //        // current_level + 1 recursion
        //        index.push(current_dim, source_level + 1, 2*source_index - 1);
        //        this->regular_rec(storage, index, current_dim, current_level + 1, level);
        //
        //        index.push(current_dim, source_level + 1, 2*source_index + 1);
        //        this->regular_rec(storage, index, current_dim, current_level + 1, level);
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
        //  void regular_rec_1d(GridStorage* storage, index_type& index, level_t current_level, level_t level)
        //  {
        //        for(level_t l = 1; l <= level-current_level + 1; l++)
        //        {
        //            if (l == level-current_level+1)
        //            {
        //              for(index_t i = 1; static_cast<int>(i) <= static_cast<int>(1<<(l-1)); i++)
        //                {
        //                    index.push(0, l, 2*i-1, true);
        //                    storage->insert(index);
        //                }
        //            }
        //            else
        //            {
        //              for(index_t i = 1; static_cast<int>(i) <= static_cast<int>(1<<(l-1)); i++)
        //                {
        //                    index.push(0, l, 2*i-1, false);
        //                    storage->insert(index);
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
        void boundaries_trapezoid_rec(GridStorage* storage, index_type& index, size_t current_dim, level_t current_level, level_t level, bool bLevelZero) {
          if (current_dim == 0) {
            boundaries_Trapezoid_rec_1d(storage, index, current_level, level, bLevelZero);
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
                this->boundaries_trapezoid_rec(storage, index, current_dim - 1, current_level, level, true);

                index.push(current_dim, 0, 1, false);
                this->boundaries_trapezoid_rec(storage, index, current_dim - 1, current_level, level, true);

                index.push(current_dim, source_level, source_index);
              }

              // d-1 recursion
              index.setLeaf(bLeafProperty);
              this->boundaries_trapezoid_rec(storage, index, current_dim - 1, current_level, level, bLevelZero);
            }

            if (current_level < level) {
              index.push(current_dim, source_level + 1, 2 * source_index - 1);
              this->boundaries_trapezoid_rec(storage, index, current_dim, current_level + 1, level, bLevelZero);

              index.push(current_dim, source_level + 1, 2 * source_index + 1);
              this->boundaries_trapezoid_rec(storage, index, current_dim, current_level + 1, level, bLevelZero);
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
        void boundaries_Trapezoid_rec_1d(GridStorage* storage, index_type& index, level_t current_level, level_t level, bool bLevelZero) {
          bool bLevelGreaterZero = !bLevelZero;

          for (level_t l = 0; l <= level - current_level + 1; l++) {
            if (l == level - current_level + 1) {
              if (l == 0) {
                index.push(0, 0, 0, false);
                storage->insert(index);
                index.push(0, 0, 1, false);
                storage->insert(index);
              } else {
                for (index_t i = 1; static_cast<int>(i) <= static_cast<int>(1 << (l - 1)); i++) {
                  index.push(0, l, 2 * i - 1, (true && bLevelGreaterZero));
                  storage->insert(index);
                }
              }
            } else {
              if (l == 0) {
                index.push(0, 0, 0, false);
                storage->insert(index);
                index.push(0, 0, 1, false);
                storage->insert(index);
              } else {
                for (index_t i = 1; static_cast<int>(i) <= static_cast<int>(1 << (l - 1)); i++) {
                  index.push(0, l, 2 * i - 1, false);
                  storage->insert(index);
                }
              }
            }
          }
        }

        /**
         * recursive construction of the spare grid with boundaries, classic level 0 approach, only for level 0 and 1
         *
         * @param storage hashmap that stores the grid points
         * @param index point's index
         * @param current_dim current working dimension
         * @param current_level current level in this construction step
         * @param level maximum level of the sparse grid
         */
        void boundaries_rec(GridStorage* storage, index_type& index, size_t current_dim, level_t current_level, level_t level) {
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
                storage->insert(index);
                index.push(0, 0, 1, bLeafProperty);
                storage->insert(index);

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
                storage->insert(index);
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
        void square_rec(GridStorage* storage, index_type& index,
                        size_t current_dim, level_t level, level_t small_level,
                        bool tail, size_t sum) {
          index_t source_index;
          level_t source_level;

          index.get(current_dim, source_level, source_index);
          bool newtail = tail || (source_level > small_level);

          // d-1 recursion
          if (source_level == 0) {
            if (current_dim == 0) {
              index.setLeaf(false);
              index.push(0, 0, 0);
              storage->insert(index);
              index.push(0, 0, 1);
              storage->insert(index);

              index.push(current_dim, source_level, source_index);
            } else {
              //we want all nodes, so we initialize the unassigned levels to start from 0
              for (size_t d = 0; d < current_dim; d++)
                index.push(d, 0, 0);

              index.push(current_dim, 0, 0);
              this->square_rec(storage, index, current_dim - 1, level, small_level, newtail, sum);


              for (size_t d = 0; d < current_dim; d++)
                index.push(d, 0, 0);

              index.push(current_dim, 0, 1);

              this->square_rec(storage, index, current_dim - 1, level, small_level, newtail, sum);

              index.push(current_dim, source_level, source_index);
            }
          } else {
            if (current_dim == 0) {
              /**
               * If a level of the node equals level, and all the others equal small_level, the node is a leaf
               * This is equivalent to saying the sum of levels equals small_level*(dim-1)+level
               * */
              if (sum == (small_level * (storage->dim() - 1) + level)) index.setLeaf(true);
              else index.setLeaf(false);

              storage->insert(index);
            } else {

              for (size_t d = 0; d < current_dim; d++)
                index.push(d, 0, 0);

              this->square_rec(storage, index, current_dim - 1,  level, small_level, newtail, sum);
            }
          }

          /**If the level of the node is smaller than small_level or we didn't have yet a level greater than small_level(!tail) and the level is smaller then level then
          we can then proceed to the next level on this dimension, otherwise we reached the maximum possible level*
          */
          if ((source_level < small_level) || ((!tail) && (source_level < level))) {
            if (source_level == 0 && source_index == 0) {
              index.push(current_dim, source_level + 1, 1);
              this->square_rec(storage, index, current_dim,  level, small_level, tail, sum + 1);
            } else {
              index.push(current_dim, source_level + 1, 2 * source_index - 1);
              this->square_rec(storage, index, current_dim,  level, small_level, tail, sum + 1);

              index.push(current_dim, source_level + 1, 2 * source_index + 1);
              this->square_rec(storage, index, current_dim,  level, small_level, tail, sum + 1);
            }
          }

        }
        /**
         * recursive construction of a super trapezoid grid with boundaries
         *
         * @param storage hashmap that stores the grid points
         * @param index point's index
         * @param current_dim current working dimension
         * @param current_level the current level of the gridpoint so far, starts from minlevel*dim
         * @param level the maximum level of the gridpoint
         * @param minlevel the level limit given by the user(tells us which fullgrids won't be present in the construction of the sparse grid)
         */
        void trunc_rec(GridStorage* storage, index_type& index,
                       size_t current_dim, level_t current_level,
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
                storage->insert(index);
                index.push(0, 0, 1, bLeafProperty);
                storage->insert(index);

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
                storage->insert(index);
              } else {
                index.setLeaf(bLeafProperty);
                trunc_rec(storage, index, current_dim - 1, current_level, level, minlevel);
              }

              index.setLeaf(bSaveLeafProperty);
            }
          }

          if (source_level < minlevel) {
            /**
             * If the source level of the node is smaller than minlevel we don't increase the variable current_level(since we started with minlevel*dim)
             * This trick makes it possible to introduce all nodes with source_level&lt;minlevel without a separate treatment
             * */
            if (source_level == 0 && source_index == 0) {
              index.push(current_dim, source_level + 1, 1);
              trunc_rec(storage, index, current_dim, current_level, level, minlevel);
            } else {
              index.push(current_dim, source_level + 1, 2 * source_index - 1);
              trunc_rec(storage, index, current_dim, current_level, level, minlevel);

              index.push(current_dim, source_level + 1, 2 * source_index + 1);
              trunc_rec(storage, index, current_dim, current_level , level, minlevel);
            }

          } else if (current_level < level) {
            /**
             * if the source_level is already >=minlevel we can proceed naturally and increase the current_level which represents the sum of levels so far
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

  }
}

#endif /* HASHGENERATOR_HPP */
