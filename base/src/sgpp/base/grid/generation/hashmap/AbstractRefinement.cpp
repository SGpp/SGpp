// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/generation/hashmap/AbstractRefinement.hpp>
#include <sgpp/base/grid/generation/refinement_strategy/RefinementDecorator.hpp>


#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace base {


/*bool refinementPairCompare(const AbstractRefinement::refinement_pair_type& element1,
                           const AbstractRefinement::refinement_pair_type& element2) {
        return element1.second > element2.second;
}
*/

size_t AbstractRefinement::getIndexOfMin(RefinementFunctor::value_type* array,
    size_t length) {
  size_t min_idx = 0;

  for (size_t i = 1; i < length; i++) {
    if (array[i] < array[min_idx])
      min_idx = i;
  }

  return min_idx;
}


void AbstractRefinement::createGridpoint1D(GridPoint& point,
    size_t d, GridStorage& storage, index_t& source_index,
    level_t& source_level) {
  point.get(d, source_level, source_index);

  if (source_level > 1) {
    if (((source_index + 1) / 2) % 2 == 1) {
      point.set(d, source_level - 1, (source_index + 1) / 2);
    } else {
      point.set(d, source_level - 1, (source_index - 1) / 2);
    }

    createGridpointSubroutine(storage, point);
    // restore values
    point.set(d, source_level, source_index);
  }
}

void AbstractRefinement::refineGridpoint1D(GridStorage& storage, size_t seq,
    size_t d) {
  this->refineGridpoint1D(storage, storage.getPoint(seq), d);
}


/*void AbstractRefinement::strategy_refine(GridStorage& storage,
        RefinementStrategy& refinement_strategy)
{
    refinement_strategy.refine(storage, this);
}*/

bool AbstractRefinement::isRefinable(GridStorage& storage, GridPoint& point) {
  GridStorage::grid_map_iterator child_iter;

  if (point.isLeaf()) return true;

  for (size_t d = 0; d < storage.getDimension(); d++) {
    index_t source_index;
    level_t source_level;
    point.get(d, source_level, source_index);

    // test existence of left child
    point.set(d, source_level + 1, 2 * source_index - 1);
    child_iter = storage.find(&point);

    // if there no more grid points --> test if we should refine the grid
    if (child_iter == storage.end()) {
      return true;
    }

    // test existance of right child
    point.set(d, source_level + 1, 2 * source_index + 1);
    child_iter = storage.find(&point);

    if (child_iter == storage.end()) {
      return true;
    }

    // reset current grid point in dimension d
    point.set(d, source_level, source_index);
  }

  return false;
}

}  // namespace base
}  // namespace sgpp

