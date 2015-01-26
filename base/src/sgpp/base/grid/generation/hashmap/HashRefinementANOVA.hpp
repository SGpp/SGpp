/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de), Dirk Pflueger (pflueged@in.tum.de)
#ifndef HASHREFINEMENT_HPP
#define HASHREFINEMENT_HPP

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/generation/functors/RefinementFunctor.hpp>

#include <sgpp/base/exception/generation_exception.hpp>

#include <vector>
#include <cmath>

#include <iostream>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * Standard free refinement class for sparse grids without boundaries
     */
    class HashRefinement {
      public:
        typedef GridStorage::index_type index_type;
        typedef index_type::index_type index_t;
        typedef index_type::level_type level_t;

        /**
         * Refines a grid according to a RefinementFunctor provided.
         * Refines all grid points in the same ANOVA component if
         * possible, and if their refinement value is larger than RefinementFunctor::start()
         * and their absolute value is larger or equal than RefinementFunctor::getRefinementThreshold()
         *
         * @param storage hashmap that stores the grid points
         * @param functor a RefinementFunctor specifying the refinement criteria
         */
        void free_refine(GridStorage* storage, RefinementFunctor* functor) {
          // only the function with the local error indicator grater than the
          // threshold will be refined
          double threshold = functor->getRefinementThreshold();

          if (storage->size() == 0) {
            throw generation_exception("storage empty");
          }

          index_type index;
          GridStorage::grid_map_iterator end_iter = storage->end();

          // start iterating over whole grid
          for (GridStorage::grid_map_iterator iter = storage->begin();
               iter != end_iter; iter++) {
            index = *(iter->first);

            GridStorage::grid_map_iterator child_iter;

            for (size_t d = 0; d < storage->dim(); d++) {
              index_t source_index;
              level_t source_level;
              index.get(d, source_level, source_index);

              // in order to remain in the same ANOVA component, we shouldn't
              // refine the constant functions (level 1)
              if (source_level == 0)
                continue;

              // test existence of the left child
              index.set(d, source_level + 1, 2 * source_index - 1);
              child_iter = storage->find(&index);

              // if there no more grid points --> test if we should refine the grid
              if (child_iter == end_iter) {
                RefinementFunctor::value_type current_value =
                  (*functor)(storage, iter->second);

                if (current_value > threshold)
                  refine_gridpoint_directional(storage, iter->second, d);

                break; // you don't have to test for the right child
              }

              // if there is a left child test the existence of the right child
              index.set(d, source_level + 1, 2 * source_index + 1);
              child_iter = storage->find(&index);

              if (child_iter == end_iter) {
                RefinementFunctor::value_type current_value =
                  (*functor)(storage, iter->second);

                if (current_value > threshold)
                  refine_gridpoint_directional(storage, iter->second, d);

                break;
              }

              // reset current grid point in dimension d
              index.set(d, source_level, source_index);
            }
          }

        }

        /**
         * Computes and returns the number of grid points, which can be refined.
         * This is the number of grid points that have at least one child missing.
         *
         * @param storage hashmap that stores the grid points
         * @return The number of grid points that can be refined
         */
        size_t getNumberOfRefinablePoints(GridStorage* storage) {
          size_t counter = 0;

          if (storage->size() == 0) {
            throw generation_exception("storage empty");
          }

          index_type index;
          GridStorage::grid_map_iterator end_iter = storage->end();

          // start iterating over whole grid
          for (GridStorage::grid_map_iterator iter = storage->begin();
               iter != end_iter; iter++) {
            index = *(iter->first);

            GridStorage::grid_map_iterator child_iter;

            // check for each grid point whether it can be refined (i.e., whether not all children exist yet)
            for (size_t d = 0; d < storage->dim(); d++) {
              index_t source_index;
              level_t source_level;
              index.get(d, source_level, source_index);

              // test existance of left child
              index.set(d, source_level + 1, 2 * source_index - 1);
              child_iter = storage->find(&index);

              // if there no more grid points --> test if we should refine the grid
              if (child_iter == end_iter) {
                counter++;
              }

              // test existance of right child
              index.set(d, source_level + 1, 2 * source_index + 1);
              child_iter = storage->find(&index);

              if (child_iter == end_iter) {
                counter++;
              }

              // reset current grid point in dimension d
              index.set(d, source_level, source_index);
            }
          }

          return counter;

        }

      protected:
        /**
         * This method refines a grid point by generating the children in every dimension
         * of the grid and all their missing ancestors by calling create_gridpoint().
         *
         * @param storage hashmap that stores the gridpoints
         * @param refine_index The index in the hashmap of the point that should be refined
         */
        void refineGridpoint(GridStorage* storage, size_t refine_index) {
          index_type index((*storage)[refine_index]);

          //Sets leaf property of index, which is refined to false
          (*storage)[refine_index]->setLeaf(false);

          // @todo (blank) Maybe it's possible to move predecessor/successor discovery into the storage concept
          for (size_t d = 0; d < storage->dim(); d++) {
            index_t source_index;
            level_t source_level;
            index.get(d, source_level, source_index);

            // generate left child, if necessary
            index.set(d, source_level + 1, 2 * source_index - 1);

            if (!storage->has_key(&index)) {
              index.setLeaf(true);
              createGridpoint(storage, index);
            }

            // generate right child, if necessary
            index.set(d, source_level + 1, 2 * source_index + 1);

            if (!storage->has_key(&index)) {
              index.setLeaf(true);
              createGridpoint(storage, index);
            }

            index.set(d, source_level, source_index);
          }
        }

        /**
         * This method refines a grid point by generating the children in the given dimension
         * of the grid and all their missing ancestors by calling create_gridpoint().
         *
         * @param storage hashmap that stores the gridpoints
         * @param refine_index The index in the hashmap of the point that should be refined
         * @param d dimension in which the grid point should be refined
         */
        void refine_gridpoint_directional(GridStorage* storage, size_t refine_index,
                                          size_t d) {
          index_type index((*storage)[refine_index]);

          //Sets leaf property of index, which is refined to false
          (*storage)[refine_index]->setLeaf(false);

          index_t source_index;
          level_t source_level;
          index.get(d, source_level, source_index);

          // generate left child, if necessary
          index.set(d, source_level + 1, 2 * source_index - 1);

          if (!storage->has_key(&index)) {
            index.setLeaf(true);
            createGridpoint(storage, index);
          }

          // generate right child, if necessary
          index.set(d, source_level + 1, 2 * source_index + 1);

          if (!storage->has_key(&index)) {
            index.setLeaf(true);
            createGridpoint(storage, index);
          }

          index.set(d, source_level, source_index);

        }

        /**
         * This method creates a new point on the grid. It checks if some parents or
         * children are needed in other dimensions.
         *
         * @param storage hashmap that stores the gridpoints
         * @param index The point that should be inserted
         */
        void createGridpoint(GridStorage* storage, index_type& index) {
          for (size_t d = 0; d < storage->dim(); d++) {
            index_t source_index;
            level_t source_level;
            index.get(d, source_level, source_index);

            if (source_level > 1) {
              // @todo (blank) Maybe it's possible to move predecessor/successor discovery into the storage concept
              if (((source_index + 1) / 2) % 2 == 1) {
                index.set(d, source_level - 1, (source_index + 1) / 2);
              } else {
                index.set(d, source_level - 1, (source_index - 1) / 2);
              }

              if (!storage->has_key(&index)) {
                // save old leaf value
                bool saveLeaf = index.isLeaf();
                index.setLeaf(false);
                createGridpoint(storage, index);
                // restore leaf value
                index.setLeaf(saveLeaf);
              } else {
                // set stored index to false
                (storage->get((storage->find(&index))->second))->setLeaf(
                  false);
              }

              // restore values
              index.set(d, source_level, source_index);
            }
          }

          storage->insert(index);
        }

        /**
         * Returns the index of the first accurance of minimal element in array.
         * Used to find which entry is to be replaced next searching the maximum ones.
         *
         * @param array array with values
         * @param length length of array
         *
         * @return index of the first occurrence of minimal element in array
         */
        size_t getIndexOfMin(RefinementFunctor::value_type* array, size_t length) {
          size_t min_idx = 0;

          for (size_t i = 1; i < length; i++) {
            if (array[i] < array[min_idx])
              min_idx = i;
          }

          return min_idx;
        }
    };

  }
}

#endif /* HASHREFINEMENT_HPP */
