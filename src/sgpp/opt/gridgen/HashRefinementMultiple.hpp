/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef HASHREFINEMENTMULTIPLE_HPP
#define HASHREFINEMENTMULTIPLE_HPP

#include "base/grid/generation/hashmap/HashRefinement.hpp"

namespace sg {
  namespace opt {
    namespace gridgen {

      /**
       * Descendant of base::HashRefinement refining without the generation of hierarchical ancestors.
       *
       * In SG++ grids (as in sg::base), every grid fulfills the "hierarchical ancestors" property,
       * e.g. every gridpoint has a direct ancestor in every dimension whose level is > 1.
       *
       * By using this refinement class, the property doesn't necessarily hold anymore.
       * The original base::HashRefinement looks at the neigbors of a grid point to be refined.
       * If a neighbor in a dimension and a direction (left/right) already exists, no new point is
       * inserted. If it doesn't exist, it is inserted and all hierarchical ancestors are generated
       * recursively.
       *
       * HashRefinementMultiple instead inserts exactly \f$2d\f$ new grid points, if a point is to
       * be refined (for Noboundary grids with \f$d\f$ being the number of dimensions).
       * This is done by looking at neighbors of higher levels (closer to the point to be refined)
       * and inserting the first neighbor that doesn't exist in the grid.
       * Additionally, no hierarchical ancestors are creating by omitting the
       * base::HashRefinement::createGridpoint() calls.
       *
       * Grids without the "hierarchical ancestors" property don't allow most standard algorithms
       * to be executed on them, therefore grids and basis functions are separated in this module
       * from those in sg::base.
       */
      class HashRefinementMultiple : public base::HashRefinement {
        public:
          /**
           * Returns the number of refineable points, e.g. the number of all points
           * (every grid point can be refined).
           *
           * @param storage   grid storage
           * @return          number of refinable points
           */
          size_t getNumberOfRefinablePoints(base::GridStorage* storage) {
            return storage->size();
          }

          /**
           * Refines a grid point in one dimension.
           * This creates exactly two new grid points.
           *
           * @param storage   grid storage
           * @param index     index of the grid point
           * @param t         dimension in which the refinement should take place
           */
          void refineGridpoint1D(base::GridStorage* storage, index_type& index, size_t t) {
            index_t source_index, child_index;
            level_t source_level, child_level;

            index.get(t, source_level, source_index);

            // don't generate the left child, if the "source" has level 0, index 0 (x = 0)
            if ((source_level > 0) || (source_index == 1)) {
              // generate left child
              child_index = source_index;
              child_level = source_level;

              while (storage->has_key(&index)) {
                child_index *= 2;
                child_level++;
                index.set(t, child_level, child_index - 1);
              }

              index.setLeaf(true);
              // instead of "createGridpoint(storage, index);"
              storage->insert(index);
              index.set(t, source_level, source_index);
            }

            // don't generate the right child, if the "source" has level 0, index 1 (x = 1)
            if ((source_level > 0) || (source_index == 0)) {
              // generate right child
              child_index = source_index;
              child_level = source_level;

              while (storage->has_key(&index)) {
                child_index *= 2;
                child_level++;
                index.set(t, child_level, child_index + 1);
              }

              index.setLeaf(true);
              // instead of "createGridpoint(storage, index);"
              storage->insert(index);
              index.set(t, source_level, source_index);
            }
          }

        protected:
          /**
           * Examine the grid points and stores the indices those that can be refined
           * and have maximal indicator values.
           * This function differences from its counterpart in base::HashRefinement
           * insofar every grid point could potentially be refined.
           *
           * @param storage           grid storage
           * @param functor           refinement criteria
           * @param refinements_num   maximal number of points to refine
           * @param max_indices       the array where the point indices should be stored
           * @param max_values        the array where the corresponding indicator values should be stored
           */
          void collectRefinablePoints(base::GridStorage* storage, base::RefinementFunctor* functor,
                                      size_t refinements_num, size_t* max_indices,
                                      base::RefinementFunctor::value_type* max_values) {
            size_t min_idx = 0;

            // max value equals min value
            base::RefinementFunctor::value_type max_value = max_values[min_idx];

            index_type index;
            base::GridStorage::grid_map_iterator end_iter = storage->end();

            // start iterating over whole grid
            for (base::GridStorage::grid_map_iterator iter = storage->begin();
                 iter != end_iter; iter++) {
              base::RefinementFunctor::value_type current_value = (*functor)(storage, iter->second);

              if (current_value > max_value) {
                // replace the minimal point in result array, find the new minimal point
                max_values[min_idx] = current_value;
                max_indices[min_idx] = iter->second;
                min_idx = getIndexOfMin(max_values, refinements_num);
                max_value = max_values[min_idx];
              }
            }
          }
      };

    }
  }
}

#endif
