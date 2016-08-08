// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/generation/refinement_strategy/ImpurityRefinement.hpp>

#include <sgpp/globaldef.hpp>
#include <sgpp/base/grid/generation/functors/ImpurityRefinementIndicator.hpp>

#include <algorithm>
#include <limits>
#include <list>
#include <vector>
#include <utility>
#include <iterator>


namespace sgpp {
namespace base {


void ImpurityRefinement::collectRefinablePoints(GridStorage& storage, 
                                                RefinementFunctor& functor, 
                                                AbstractRefinement::refinement_container_type& collection) {
  size_t refinements_num = functor.getRefinementsNum();
  GridStorage::grid_map_iterator end_iter = storage.end();

  // start iterating over whole grid
  for (GridStorage::grid_map_iterator iter = storage.begin(); iter != end_iter; iter++) {
    AbstractRefinement::refinement_list_type current_value_list = getIndicator(storage, iter, functor);
    addElementToCollection(iter, current_value_list, refinements_num, collection);
  }
}

AbstractRefinement::refinement_list_type ImpurityRefinement::getIndicator(
  GridStorage& storage,
  const GridStorage::grid_map_iterator& iter,
  const RefinementFunctor& functor) const {
  AbstractRefinement::refinement_list_type list;

  // this refinement algorithm uses the impurity refinement indicator.
  // dynamic casting is used to maintain the signature of the algorithm,
  // but still be able to use the
  // impurity refinement indicator with it.
  const ImpurityRefinementIndicator& impurityIndicator =
    dynamic_cast<const ImpurityRefinementIndicator&>(functor);
  refinement_key_type* key;

  GridPoint& point = *(iter->first);
  GridStorage::grid_map_iterator child_iter;
  GridStorage::grid_map_iterator end_iter = storage.end();

  bool isRefinablePoint = false;

  if (point.isLeaf()) {
    isRefinablePoint = true;
  }

  for (size_t d = 0; d < storage.getDimension(); d++) {
    index_t source_index;
    level_t source_level;
    point.get(d, source_level, source_index);

    // test existence of left child
    point.set(d, source_level + 1, 2 * source_index - 1);
    child_iter = storage.find(&point);

    // if there no more grid points --> test if we should refine the grid
    if (child_iter == storage.end()) {
      isRefinablePoint = true;
      point.set(d, source_level, source_index);
      break;
    }

    // test existance of right child
    point.set(d, source_level + 1, 2 * source_index + 1);
    child_iter = storage.find(&point);

    if (child_iter == storage.end()) {
      isRefinablePoint = true;
      point.set(d, source_level, source_index);
      break;
    }

    // reset current grid point in dimension d
    point.set(d, source_level, source_index);
  }

  if (isRefinablePoint) {
    
    double impurity = impurityIndicator(point);

    if (impurity > iThreshold_) {
      size_t d = 0;
      key = new refinement_key_type(point, storage.getSequenceNumber(point), d);        // d actually not required now
      list.emplace_front(std::shared_ptr<AbstractRefinement::refinement_key_type>(key),
                         impurity); 
    }
  }

  return list;
}


void ImpurityRefinement::addElementToCollection(
  const GridStorage::grid_map_iterator& iter,
  AbstractRefinement::refinement_list_type current_value_list,
  size_t refinements_num,
  AbstractRefinement::refinement_container_type& collection) {
  for (AbstractRefinement::refinement_list_type::iterator it =
         current_value_list.begin(); it != current_value_list.end(); it++) {
    collection.push_back(*it);
    std::push_heap(collection.begin(), collection.end(), AbstractRefinement::compare_pairs);

    if (collection.size() > refinements_num) {
      // remove the top (smallest) element
      std::pop_heap(collection.begin(), collection.end(),
                    AbstractRefinement::compare_pairs);
      collection.pop_back();
    }

  }
}


void ImpurityRefinement::refineGridpointsCollection(
  GridStorage& storage, RefinementFunctor& functor,
  AbstractRefinement::refinement_container_type& collection) {

  // now refine all grid points which satisfy the refinement criteria

  //double threshold = functor.getRefinementThreshold();
  refinement_key_type* key;

  ImpurityRefinementIndicator& impurityIndicator =
    dynamic_cast<ImpurityRefinementIndicator&>(functor);

  for (AbstractRefinement::refinement_pair_type& pair : collection) {
    key = dynamic_cast<refinement_key_type*>(pair.first.get());

    GridPoint& point = key->getPoint();

    //std::cout << "refine point: " << storage.getSequenceNumber(point) << std::endl;

    storage[storage.getSequenceNumber(point)].setLeaf(false);
    //point.setLeaf(false);

    for (size_t d = 0; d < storage.getDimension(); d++) {
      index_t source_index;
      level_t source_level;
      point.get(d, source_level, source_index);
      // generate left child, if necessary
      point.set(d, source_level + 1, 2 * source_index - 1);

      if (!storage.isContaining(point)) {
        point.setLeaf(true);
        createGridpointMod(storage, point, impurityIndicator);
      }

      // generate right child, if necessary
      point.set(d, source_level + 1, 2 * source_index + 1);
 
      if (!storage.isContaining(point)) {
        point.setLeaf(true);
        createGridpointMod(storage, point, impurityIndicator);
      }
 
      point.set(d, source_level, source_index);
    }
    
  }

  collection.empty();
}


void ImpurityRefinement::createGridpointMod(
  GridStorage& storage, GridPoint& point, ImpurityRefinementIndicator& indicator) {
  index_t source_index;
  level_t source_level;

  for (size_t d = 0; d < storage.getDimension(); d++) {
    createGridpoint1DMod(point, d, storage, source_index, source_level, indicator);
  }

  storage.insert(point);
  indicator.update(point);
}


void ImpurityRefinement::createGridpoint1DMod(GridPoint& point, size_t d, GridStorage& storage,
                                                     index_t& source_index, level_t& source_level,
                                                     ImpurityRefinementIndicator& indicator) {
  point.get(d, source_level, source_index);

  if (source_level > 1) {
    if (((source_index + 1) / 2) % 2 == 1) {
      point.set(d, source_level - 1, (source_index + 1) / 2);
    } else {
      point.set(d, source_level - 1, (source_index - 1) / 2);
    }

    // grid point subroutine
    if (!storage.isContaining(point)) {
      // save old leaf value
      bool saveLeaf = point.isLeaf();
      point.setLeaf(false);
      createGridpointMod(storage, point, indicator);
      // restore leaf value
      point.setLeaf(saveLeaf);
    } else {
      // set stored index to false
      storage.getPoint((storage.find(&point))->second).setLeaf(false);
    }
    // restore values
    point.set(d, source_level, source_index);
  }
}


void ImpurityRefinement::free_refine(GridStorage& storage, ImpurityRefinementIndicator& functor) {
  if (storage.getSize() == 0) {
    throw generation_exception("storage empty");
  }

  // initialization
  AbstractRefinement::refinement_container_type collection;

  collectRefinablePoints(storage, functor, collection);
  // now refine all grid points which satisfy the refinement criteria
  refineGridpointsCollection(storage, functor, collection);
  collection.clear();
}


}  // namespace base
}  // namespace sgpp
