// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/generation/refinement_strategy/SubspaceRefinement.hpp>
#include <sgpp/base/grid/generation/hashmap/AbstractRefinement.hpp>
#include <sgpp/base/grid/common/IndexInSubspaceGenerator.hpp>

#include <sgpp/globaldef.hpp>

#include <algorithm>
#include <iostream>
#include <vector>

namespace sgpp {
namespace base {

void SubspaceRefinement::addElementToCollection(
  const GridStorage::grid_map_iterator& iter,
  AbstractRefinement::refinement_list_type current_value_list,
  size_t refinements_num,
  AbstractRefinement::refinement_container_type& collection) {
  for (AbstractRefinement::refinement_pair_type key_value :
       current_value_list) {
    if (key_value.second == 0) continue;

    auto predicat = [key_value](const AbstractRefinement::refinement_pair_type &
    element) {
      bool result = true;

      for (size_t d = 0; d < key_value.first->getPoint().getDimension(); d++) {
        result = result
                 && (key_value.first->getPoint().getLevel(d) ==
                     element.first->getPoint().getLevel(d));
      }

      return result;
    };

    AbstractRefinement::refinement_container_type::iterator iter = std::find_if(
          collection.begin(),
          collection.end(), predicat);

    if (iter != collection.end()) {
      // subspace with this level is already in collection,
      // hence increase the value
      iter->second += key_value.second;
    } else {
      // subspace is not yet in the collection, hence add it
      collection.push_back(key_value);
    }
  }
}


void SubspaceRefinement::free_refine(GridStorage& storage,
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



void SubspaceRefinement::collectRefinablePoints(GridStorage& storage,
    RefinementFunctor& functor,
    AbstractRefinement::refinement_container_type& collection) {
  if (functor.getRefinementsNum() == 0) return;

  size_t refinements_num = functor.getRefinementsNum();

  GridPoint point;
  GridStorage::grid_map_iterator end_iter = storage.end();


  // start iterating over whole grid
  for (GridStorage::grid_map_iterator iter = storage.begin();
       iter != end_iter; iter++) {
    point = *(iter->first);

    // std::cout <<"grid point " << iter->second << std::endl;

    GridStorage::grid_map_iterator child_iter;

    // check for each grid point whether it can be refined
    // (i.e., whether not all kids exist yet)
    // if yes, check whether it belongs to the refinements_num largest ones
    for (size_t d = 0; d < storage.getDimension(); d++) {
      index_t source_index;
      level_t source_level;
      point.get(d, source_level, source_index);

      // std::cout << "processing dimension " << d << std::endl;

      // test existence of left child
      point.set(d, source_level + 1, 2 * source_index - 1);
      child_iter = storage.find(&point);

      // if there no more grid points --> test if we should refine the grid
      if (child_iter == end_iter) {
        AbstractRefinement::refinement_list_type current_value_list =
          getIndicator(storage, iter,
                       functor);

        // std::cout << "size of the list " << current_value_list.empty() <<
        // std::endl;

        addElementToCollection(iter, current_value_list,
                               refinements_num, collection);
        break;
      }

      // test existence of right child
      point.set(d, source_level + 1, 2 * source_index + 1);
      child_iter = storage.find(&point);

      if (child_iter == end_iter) {
        AbstractRefinement::refinement_list_type current_value_list =
          getIndicator(storage, iter,
                       functor);

        addElementToCollection(iter, current_value_list,
                               refinements_num, collection);
        break;
      }

      // reset current grid point in dimension d
      point.set(d, source_level, source_index);
    }
  }




  // nth_element makes sure that the first refinements_num elements in the
  // vector are larger then the rest
  std::nth_element(collection.begin(),
                   collection.begin() + refinements_num,
                   collection.end(), AbstractRefinement::compare_pairs);

  if (collection.size() >= refinements_num) {
    // clear the collection and populate it only with those elements
    // that will be refined
    for (size_t diff = collection.size() - refinements_num; diff > 0; diff--) {
      collection.pop_back();
    }
  }
}


void SubspaceRefinement::refineGridpointsCollection(GridStorage& storage,
    RefinementFunctor& functor,
    AbstractRefinement::refinement_container_type& collection) {

  HashGridPoint grid_index(storage.getDimension());

  // refine all points of the subspace in all dimensions
  for (AbstractRefinement::refinement_pair_type& pair : collection) {
    const std::vector<level_t> level_vector = pair.first->getLevelVector();
    IndexInSubspaceGenerator subspace(level_vector);

    for (IndexInSubspaceGenerator::iterator index_it = subspace.begin();
         index_it != subspace.end(); index_it++) {
      for (size_t d = 0; d < storage.getDimension(); d++) {
        grid_index.set(d, level_vector[d], (*index_it)[d]);
      }

      const size_t seq = storage.getSequenceNumber(grid_index);

      if (seq < storage.getSize()) {
        refineGridpoint(storage, seq);
      }
    }
  }
}

}  // namespace base
}  // namespace sgpp
