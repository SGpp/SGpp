// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/generation/hashmap/HashRefinementInteraction.hpp>

#include <algorithm>
#include <memory>
#include <vector>

namespace sgpp {
namespace base {

HashRefinementInteraction::HashRefinementInteraction(
    std::unordered_set<std::vector<bool>> interactions)
    : HashRefinement(), interactions(interactions) {}

void HashRefinementInteraction::createGridpoint(GridStorage& storage, GridPoint& index) {
  index_t source_index;
  level_t source_level;

  // Get the currently used dimensions of the grid point.
  auto coordsIndex = std::vector<bool>(storage.getDimension());
  for (size_t i = 0; i < storage.getDimension(); ++i) {
    coordsIndex[i] = index.getStandardCoordinate(i) != 0.5;
  }

  // Check first whether index obeys our limitations.
  if (interactions.find(coordsIndex) == interactions.end()) {
    return;
  }

  // We now check whether the index with additional dimension d is legal.
  for (size_t d = 0; d < storage.getDimension(); d++) {
    auto curCoords = coordsIndex;
    curCoords[d] = true;
    if (interactions.find(curCoords) != interactions.end()) {
      createGridpoint1D(index, d, storage, source_index, source_level);
    }
  }

  storage.insert(index);
}

void HashRefinementInteraction::collectRefinablePoints(GridStorage& storage,
    RefinementFunctor& functor,
    AbstractRefinement::refinement_container_type& collection) {
  size_t refinements_num = functor.getRefinementsNum();

  // max value equals min value

  GridPoint point;
  GridStorage::grid_map_iterator end_iter = storage.end();

  // start iterating over whole grid
  for (GridStorage::grid_map_iterator iter = storage.begin(); iter != end_iter;
       iter++) {
    point = *(iter->first);

    // get the boolean interactions of current point
    auto coordsIndex = std::vector<bool>(storage.getDimension());
    for (size_t i = 0; i < storage.getDimension(); ++i) {
      coordsIndex[i] = point.getStandardCoordinate(i) != 0.5;
    }

    GridStorage::grid_map_iterator child_iter;

    // check for each grid point whether it can be refined
    // (i.e., whether not all kids exist yet)
    // if yes, check whether it belongs to the refinements_num largest ones
    for (size_t d = 0; d < storage.getDimension(); d++) {
      auto tmpCoord = coordsIndex[d];
      coordsIndex[d] = 1;
      // point isnt included in interactionterms
      if (interactions.find(coordsIndex) == interactions.end()) {
        coordsIndex[d] = tmpCoord;
        break;
      }
      index_t source_index;
      level_t source_level;
      point.get(d, source_level, source_index);

      // test existence of left child
      point.set(d, source_level + 1, 2 * source_index - 1);
      child_iter = storage.find(&point);

      // if there no more grid points --> test if we should refine the grid
      if (child_iter == end_iter) {
        AbstractRefinement::refinement_list_type current_value_list =
          getIndicator(storage, iter, functor);
        addElementToCollection(iter, current_value_list, refinements_num,
                               collection);
        break;
      }

      // test existence of right child
      point.set(d, source_level + 1, 2 * source_index + 1);
      child_iter = storage.find(&point);

      if (child_iter == end_iter) {
        AbstractRefinement::refinement_list_type current_value_list =
          getIndicator(storage, iter, functor);
        addElementToCollection(iter, current_value_list, refinements_num,
                               collection);
        break;
      }

      // reset current grid point in dimension d
      point.set(d, source_level, source_index);
    }
  }
}
}  // namespace base
}  // namespace sgpp
