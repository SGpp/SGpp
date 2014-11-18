/*
 * RefinementANOVAStrategy.cpp
 *
 *  Created on: Mar 6, 2012
 *      Author: khakhutv_local
 */

//#include "base/grid/generation/hashmap/HashRefinementAbstract.hpp"
#include "base/grid/generation/refinement_strategy/RefinementANOVAStrategy.hpp"
#include "base/grid/generation/functors/RefinementFunctor.hpp"
#include <list>

using namespace std;

namespace sg {
  namespace base {


    IndexDimension RefinementANOVAStrategy::createIndexDimensionItem(HashRefinementAbstract::index_type* index, size_t d) {
      IndexDimension item;
      HashRefinementAbstract::index_type i = new HashRefinementAbstract::index_type(index);
      item.d = d;
      item.index = &i;
      return item;
    }

    void RefinementANOVAStrategy::refine(GridStorage* storage,
                                         HashRefinementAbstract* hash_refinement) {
      // only the function with the local error indicator grater than the
      // threshold will be refined
      double threshold = get_refinement_functor()->getRefinementThreshold()
                         * get_refinement_functor()->getTotalRefinementValue(storage);
      vector<HashRefinementAbstract::index_type*> index_vector;
      vector<size_t> dim_vector;

      if (storage->size() == 0) {
        throw generation_exception("storage empty");
      }

      HashRefinementAbstract::index_type index, index_child;
      GridStorage::grid_map_iterator end_iter = storage->end();

      // start iterating over whole grid
      for (GridStorage::grid_map_iterator iter = storage->begin();
           iter != end_iter; iter++) {
        index_child = index = *(iter->first);

        GridStorage::grid_map_iterator child_iter;

        RefinementFunctor::value_type current_value = fabs(
              (*get_refinement_functor())(storage, iter->second));

        if (current_value > threshold) {

          for (size_t d = 0; d < storage->dim(); d++) {
            HashRefinementAbstract::index_t source_index;
            HashRefinementAbstract::level_t source_level;
            index.get(d, source_level, source_index);

            // in order to remain in the same ANOVA component, we shouldn't
            // refine the constant functions (level 1)
            if (source_level <= 1)
              continue;

            // test existence of the left child
            index_child.set(d, source_level + 1, 2 * source_index - 1);
            child_iter = storage->find(&index_child);

            // if there no more grid points --> test if we should refine the grid
            if (child_iter == end_iter) {
              index_vector.push_back(new HashRefinementAbstract::index_type(index.getPointer()));
              dim_vector.push_back(d);
            } else {
              // if there is a left child test the existence of the right child
              index_child.set(d, source_level + 1, 2 * source_index + 1);
              child_iter = storage->find(&index_child);

              if (child_iter == end_iter) {
                index_vector.push_back(new HashRefinementAbstract::index_type(index.getPointer()));
                dim_vector.push_back(d);
              }
            }

            // reset current grid point in dimension d
            //index.set(d, source_level, source_index);
          }
        }
      }

      for (size_t i = 0; i < index_vector.size(); i++) {
        hash_refinement->refine_gridpoint_1d(storage, *index_vector[i], dim_vector[i]);
      }
    }

  } /* namespace base */
} /* namespace sg */
