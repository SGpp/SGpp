// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/generation/hashmap/HashRefinement.hpp>
#include <sgpp/base/exception/generation_exception.hpp>

#include <sgpp/globaldef.hpp>

#include <vector>
#include <algorithm>
#include <memory>


namespace sgpp {
namespace base {

void HashRefinement::addElementToCollection(
  const GridStorage::grid_map_iterator& iter,
  AbstractRefinement::refinement_list_type current_value_list,
  size_t refinements_num,
  AbstractRefinement::refinement_container_type& collection) {
  for (AbstractRefinement::refinement_list_type::iterator it =
         current_value_list.begin();
       it != current_value_list.end(); it++) {
    collection.push_back(*it);
    std::push_heap(collection.begin(), collection.end(),
                   AbstractRefinement::compare_pairs);


    if (collection.size() > refinements_num) {
      // remove the top (smallest) element
      std::pop_heap(collection.begin(), collection.end(),
                    AbstractRefinement::compare_pairs);
      collection.pop_back();
    }
  }
}


void HashRefinement::collectRefinablePoints(GridStorage& storage,
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

    GridStorage::grid_map_iterator child_iter;

    // check for each grid point whether it can be refined
    // (i.e., whether not all kids exist yet)
    // if yes, check whether it belongs to the refinements_num largest ones
    for (size_t d = 0; d < storage.getDimension(); d++) {
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

AbstractRefinement::refinement_list_type HashRefinement::getIndicator(
  GridStorage& storage,
  const GridStorage::grid_map_iterator& iter,
  const RefinementFunctor& functor) const {
  AbstractRefinement::refinement_list_type list;
  list.emplace_front(std::make_shared<AbstractRefinement::refinement_key_type>(*
                     (iter->first), iter->second),
                     functor(storage, iter->second));
  return list;
}



void HashRefinement::refineGridpointsCollection(GridStorage& storage,
    RefinementFunctor& functor,
    AbstractRefinement::refinement_container_type& collection) {

  double threshold = functor.getRefinementThreshold();

  for (AbstractRefinement::refinement_pair_type& pair : collection) {
    if (pair.second >= threshold) {
      refineGridpoint(storage, pair.first->getSeq());
    }
  }
}

void HashRefinement::free_refine(GridStorage& storage,
                                 RefinementFunctor& functor,
                                 std::vector<size_t>* addedPoints) {
  if (storage.getSize() == 0) {
    throw generation_exception("storage empty");
  }
  /**
   * Assumption: during the refinement process the only change made
   * to the storage is the following:
   * New (if any) gridpoints are appended (to the end) of the storage
   */
  size_t sizeBeforeRefine = storage.getSize();

  AbstractRefinement::refinement_container_type collection;
  collectRefinablePoints(storage, functor, collection);
  // now refine all grid points which satisfy the refinement criteria
  refineGridpointsCollection(storage, functor, collection);

  if (addedPoints != nullptr) {
    for (size_t i = sizeBeforeRefine; i < storage.getSize(); i++) {
      addedPoints->push_back(i);
    }
  }
}

size_t HashRefinement::getNumberOfRefinablePoints(GridStorage& storage) {
  size_t counter = 0;

  if (storage.getSize() == 0) {
    throw generation_exception("storage empty");
  }

  GridPoint point;
  GridStorage::grid_map_iterator end_iter = storage.end();

  // start iterating over whole grid
  for (GridStorage::grid_map_iterator iter = storage.begin(); iter != end_iter;
       iter++) {
    point = *(iter->first);

    GridStorage::grid_map_iterator child_iter;

    // check for each grid point whether it can be refined
    // (i.e., whether not all children exist yet)
    for (size_t d = 0; d < storage.getDimension(); d++) {
      index_t source_index;
      level_t source_level;
      point.get(d, source_level, source_index);

      // test existence of the left child
      point.set(d, source_level + 1, 2 * source_index - 1);
      child_iter = storage.find(&point);

      // if there no more grid points --> test if we should refine the grid
      if (child_iter == end_iter) {
        counter++;
        break;
      }

      // test existence of the right child
      point.set(d, source_level + 1, 2 * source_index + 1);
      child_iter = storage.find(&point);

      if (child_iter == end_iter) {
        counter++;
        break;
      }

      // reset current grid point in dimension d
      point.set(d, source_level, source_index);
    }
  }

  return counter;
}

void HashRefinement::refineGridpoint1D(GridStorage& storage, GridPoint& point, size_t d) {
  index_t source_index;
  level_t source_level;
  point.get(d, source_level, source_index);
  // generate left child, if necessary
  point.set(d, source_level + 1, 2 * source_index - 1);

  if (!storage.isContaining(point)) {
    point.setLeaf(true);
    createGridpoint(storage, point);
  }

  // generate right child, if necessary
  point.set(d, source_level + 1, 2 * source_index + 1);

  if (!storage.isContaining(point)) {
    point.setLeaf(true);
    createGridpoint(storage, point);
  }

  point.set(d, source_level, source_index);
}

void HashRefinement::refineGridpoint1D(GridStorage& storage, size_t seq,
                                       size_t d) {
  this->refineGridpoint1D(storage, storage.getPoint(seq), d);
}

void HashRefinement::refineGridpoint(GridStorage& storage,
                                     size_t refine_index) {
  GridPoint point(storage[refine_index]);
  // Sets leaf property of index, which is refined to false
  storage[refine_index].setLeaf(false);

  for (size_t d = 0; d < storage.getDimension(); d++) {
    refineGridpoint1D(storage, point, d);
  }
}

void HashRefinement::createGridpoint(GridStorage& storage, GridPoint& point) {
  index_t source_index;
  level_t source_level;

  for (size_t d = 0; d < storage.getDimension(); d++) {
    createGridpoint1D(point, d, storage, source_index, source_level);
  }

  storage.insert(point);
}

}  // namespace base
}  // namespace sgpp
