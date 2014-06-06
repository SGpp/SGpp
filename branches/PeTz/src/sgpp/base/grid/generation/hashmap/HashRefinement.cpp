/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de), Dirk Pflueger (pflueged@in.tum.de)

#include "base/grid/generation/hashmap/HashRefinement.hpp"


#include "base/exception/generation_exception.hpp"

#include <vector>
#include <cmath>
#include <iostream>


namespace sg {
  namespace base {


    void HashRefinement::collectRefinablePoints(GridStorage* storage, RefinementFunctor* functor, size_t refinements_num, size_t* max_indices, RefinementFunctor::value_type* max_values) {
      size_t min_idx = 0;

      // max value equals min value
      RefinementFunctor::value_type max_value = max_values[min_idx];
      //size_t max_index = max_indices[min_idx];

      index_type index;
      GridStorage::grid_map_iterator end_iter = storage->end();

      // start iterating over whole grid
      for (GridStorage::grid_map_iterator iter = storage->begin(); iter != end_iter; iter++) {
        index = *(iter->first);

        GridStorage::grid_map_iterator child_iter;

        // @todo (blank) Maybe it's possible to move predecessor/successor discovery into the storage concept
        // check for each grid point whether it can be refined (i.e., whether not all kids exist yet)
        // if yes, check whether it belongs to the refinements_num largest ones
        for (size_t d = 0; d < storage->dim(); d++) {
          index_t source_index;
          level_t source_level;
          index.get(d, source_level, source_index);

          // test existence of left child
          index.set(d, source_level + 1, 2 * source_index - 1);
          child_iter = storage->find(&index);

          // if there no more grid points --> test if we should refine the grid
          if (child_iter == end_iter) {
            RefinementFunctor::value_type current_value = (*functor)(storage, iter->second);

            if (current_value > max_value) {
              // replace the minimal point in result array, find the new  minimal point
              max_values[min_idx] = current_value;
              max_indices[min_idx] = iter->second;
              min_idx = getIndexOfMin(max_values, refinements_num);
              max_value = max_values[min_idx];
              break;
            }
          }

          // test existance of right child
          index.set(d, source_level + 1, 2 * source_index + 1);
          child_iter = storage->find(&index);

          if (child_iter == end_iter) {
            RefinementFunctor::value_type current_value = (*functor)(storage, iter->second);

            if (current_value > max_value) {
              // replace the minimal point in result array, find the new minimal point
              max_values[min_idx] = current_value;
              max_indices[min_idx] = iter->second;
              min_idx = getIndexOfMin(max_values, refinements_num);
              max_value = max_values[min_idx];
              break;
            }
          }

          // reset current grid point in dimension d
          index.set(d, source_level, source_index);
        }
      }
    }

    void HashRefinement::refineGridpointsCollection(GridStorage* storage, RefinementFunctor* functor, size_t refinements_num, size_t* max_indices, RefinementFunctor::value_type* max_values) {
      RefinementFunctor::value_type max_value;
      size_t max_index;
      // now refine all grid points which satisfy the refinement criteria
      double threshold = functor->getRefinementThreshold();

      for (size_t i = 0; i < refinements_num; i++) {
        max_value = max_values[i];
        max_index = max_indices[i];

        if (max_value > functor->start() && fabs(max_value) >= threshold) {
          refineGridpoint(storage, max_index);
        }
      }
    }

    void HashRefinement::free_refine(GridStorage* storage, RefinementFunctor* functor) {
      if (storage->size() == 0) {
        throw generation_exception("storage empty");
      }

      // the functor->getRefinementsNum() largest grid points should be refined.
      // gather them in an array max_values
      size_t refinements_num = functor->getRefinementsNum();
      // values
      RefinementFunctor::value_type* max_values = new RefinementFunctor::value_type[refinements_num];
      // indices
      size_t* max_indices = new size_t [refinements_num];

      // initialization
      for (size_t i = 0; i < refinements_num; i++) {
        max_values[i] = functor->start();
        max_indices[i] = 0;
      }

      collectRefinablePoints(storage, functor, refinements_num, max_indices, max_values);
      // now refine all grid points which satisfy the refinement criteria
      refineGridpointsCollection(storage, functor, refinements_num, max_indices, max_values);
      delete [] max_values;
      delete [] max_indices;

    }


    size_t HashRefinement::getNumberOfRefinablePoints(GridStorage* storage) {
      size_t counter = 0;

      if (storage->size() == 0) {
        throw generation_exception("storage empty");
      }

      index_type index;
      GridStorage::grid_map_iterator end_iter = storage->end();

      // start iterating over whole grid
      for (GridStorage::grid_map_iterator iter = storage->begin(); iter != end_iter; iter++) {
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
            break;
          }

          // test existance of right child
          index.set(d, source_level + 1, 2 * source_index + 1);
          child_iter = storage->find(&index);

          if (child_iter == end_iter) {
            counter++;
            break;
          }

          // reset current grid point in dimension d
          index.set(d, source_level, source_index);
        }
      }

      return counter;

    }


    void HashRefinement::refineGridpoint1D(GridStorage* storage, index_type& index, size_t d) {
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


    void HashRefinement::refineGridpoint(GridStorage* storage, size_t refine_index) {
      index_type index((*storage)[refine_index]);
      //Sets leaf property of index, which is refined to false
      (*storage)[refine_index]->setLeaf(false);

      // @todo (blank) Maybe it's possible to move predecessor/successor discovery into the storage concept
      for (size_t d = 0; d < storage->dim(); d++) {
        refineGridpoint1D(storage, index, d);
      }
    }




    void HashRefinement::createGridpoint(GridStorage* storage, index_type& index) {
      index_t source_index;
      level_t source_level;

      for (size_t d = 0; d < storage->dim(); d++) {
        createGridpoint1D(index, d, storage, source_index, source_level);
      }

      storage->insert(index);
    }

  }
}
