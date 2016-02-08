// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/generation/hashmap/HashRefinement.hpp>
#include <sgpp/base/exception/generation_exception.hpp>

#include <sgpp/globaldef.hpp>

#include <algorithm>
#include <memory>


namespace SGPP {
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


void HashRefinement::collectRefinablePoints(GridStorage* storage,
    RefinementFunctor* functor,
    AbstractRefinement::refinement_container_type& collection) {
  size_t refinements_num = functor->getRefinementsNum();

  // max value equals min value

  index_type index;
  GridStorage::grid_map_iterator end_iter = storage->end();

  // start iterating over whole grid
  for (GridStorage::grid_map_iterator iter = storage->begin(); iter != end_iter;
       iter++) {
    index = *(iter->first);

    GridStorage::grid_map_iterator child_iter;

    // check for each grid point whether it can be refined
    // (i.e., whether not all kids exist yet)
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
        AbstractRefinement::refinement_list_type current_value_list =
          getIndicator(storage, iter, functor);
        addElementToCollection(iter, current_value_list, refinements_num,
                               collection);
        break;
      }

      // test existence of right child
      index.set(d, source_level + 1, 2 * source_index + 1);
      child_iter = storage->find(&index);

      if (child_iter == end_iter) {
        AbstractRefinement::refinement_list_type current_value_list =
          getIndicator(storage, iter, functor);
        addElementToCollection(iter, current_value_list, refinements_num,
                               collection);
        break;
      }

      // reset current grid point in dimension d
      index.set(d, source_level, source_index);
    }
  }
}

AbstractRefinement::refinement_list_type HashRefinement::getIndicator(
  GridStorage* storage,
  const GridStorage::grid_map_iterator& iter,
  const RefinementFunctor* functor) const {
  AbstractRefinement::refinement_list_type list;
  list.emplace_front(std::make_shared<AbstractRefinement::refinement_key_type>(*
                     (iter->first), iter->second),
                     (*functor)(storage, iter->second));
  return list;
}



void HashRefinement::refineGridpointsCollection(GridStorage* storage,
    RefinementFunctor* functor,
    AbstractRefinement::refinement_container_type& collection) {

  float_t threshold = functor->getRefinementThreshold();

  for (AbstractRefinement::refinement_pair_type& pair : collection) {
    if (pair.second >= threshold) {
      refineGridpoint(storage, pair.first->getSeq());
    }
  }
}

void HashRefinement::free_refine(GridStorage* storage,
                                 RefinementFunctor* functor) {
  if (storage->size() == 0) {
    throw generation_exception("storage empty");
  }

  AbstractRefinement::refinement_container_type collection;
  collectRefinablePoints(storage, functor, collection);
  // now refine all grid points which satisfy the refinement criteria
  refineGridpointsCollection(storage, functor, collection);
}

size_t HashRefinement::getNumberOfRefinablePoints(GridStorage* storage) {
  size_t counter = 0;

  if (storage->size() == 0) {
    throw generation_exception("storage empty");
  }

  index_type index;
  GridStorage::grid_map_iterator end_iter = storage->end();

  // start iterating over whole grid
  for (GridStorage::grid_map_iterator iter = storage->begin(); iter != end_iter;
       iter++) {
    index = *(iter->first);

    GridStorage::grid_map_iterator child_iter;

    // check for each grid point whether it can be refined
    // (i.e., whether not all children exist yet)
    for (size_t d = 0; d < storage->dim(); d++) {
      index_t source_index;
      level_t source_level;
      index.get(d, source_level, source_index);

      // test existence of the left child
      index.set(d, source_level + 1, 2 * source_index - 1);
      child_iter = storage->find(&index);

      // if there no more grid points --> test if we should refine the grid
      if (child_iter == end_iter) {
        counter++;
        break;
      }

      // test existence of the right child
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

void HashRefinement::refineGridpoint1D(GridStorage* storage, index_type& index,
                                       size_t d) {
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

void HashRefinement::refineGridpoint1D(GridStorage* storage, size_t seq,
                                       size_t d) {
  this->refineGridpoint1D(storage, *(storage->get(seq)), d);
}

void HashRefinement::refineGridpoint(GridStorage* storage,
                                     size_t refine_index) {
  index_type index(*(*storage)[refine_index]);
  // Sets leaf property of index, which is refined to false
  (*storage)[refine_index]->setLeaf(false);

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

}  // namespace base
}  // namespace SGPP
