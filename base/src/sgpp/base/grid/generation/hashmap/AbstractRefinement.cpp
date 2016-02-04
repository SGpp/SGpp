// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/generation/hashmap/AbstractRefinement.hpp>
#include <sgpp/base/grid/generation/refinement_strategy/RefinementDecorator.hpp>


#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {


    /*bool refinementPairCompare(const AbstractRefinement::refinement_pair_type& element1,
                               const AbstractRefinement::refinement_pair_type& element2) {
            return element1.second > element2.second;
    }
*/

    size_t AbstractRefinement::getIndexOfMin(RefinementFunctor::value_type* array, size_t length) {
      size_t min_idx = 0;

      for (size_t i = 1; i < length; i++) {
        if (array[i] < array[min_idx])
          min_idx = i;
      }

      return min_idx;
    }


    void AbstractRefinement::createGridpoint1D(index_type& index,
        size_t d, GridStorage* storage, index_t& source_index, level_t& source_level) {
      index.get(d, source_level, source_index);

      if (source_level > 1) {
        if (((source_index + 1) / 2) % 2 == 1) {
          index.set(d, source_level - 1, (source_index + 1) / 2);
        } else {
          index.set(d, source_level - 1, (source_index - 1) / 2);
        }

        createGridpointSubroutine(storage, index);
        // restore values
        index.set(d, source_level, source_index);
      }
    }

    void AbstractRefinement::refineGridpoint1D(GridStorage* storage, size_t seq, size_t d) {
      this->refineGridpoint1D(storage, *(storage->get(seq)), d);
    }


    /*void AbstractRefinement::strategy_refine(GridStorage* storage,
            RefinementStrategy& refinement_strategy)
    {
        refinement_strategy.refine(storage, this);
    }*/

    bool AbstractRefinement::isRefinable(GridStorage* storage, index_type& index) {
      GridStorage::grid_map_iterator child_iter;

      if (index.isLeaf()) return true;

      for (size_t d = 0; d < storage->dim(); d++) {
        index_t source_index;
        level_t source_level;
        index.get(d, source_level, source_index);

        // test existence of left child
        index.set(d, source_level + 1, 2 * source_index - 1);
        child_iter = storage->find(&index);

        // if there no more grid points --> test if we should refine the grid
        if (child_iter == storage->end()) {
          return true;
        }

        // test existance of right child
        index.set(d, source_level + 1, 2 * source_index + 1);
        child_iter = storage->find(&index);

        if (child_iter == storage->end()) {
          return true;
        }

        // reset current grid point in dimension d
        index.set(d, source_level, source_index);
      }

      return false;
    }

  }
}

