// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/generation/refinement_strategy/ForwardSelectorRefinement.hpp>

#include <sgpp/globaldef.hpp>
#include <sgpp/base/grid/generation/functors/ForwardSelectorRefinementIndicator.hpp>

#include <algorithm>
#include <limits>
#include <list>
#include <vector>
#include <utility>
#include <iterator>


namespace sgpp {
namespace base {


void ForwardSelectorRefinement::collectRefinablePoints(GridStorage& storage, 
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

AbstractRefinement::refinement_list_type ForwardSelectorRefinement::getIndicator(
  GridStorage& storage,
  const GridStorage::grid_map_iterator& iter,
  const RefinementFunctor& functor) const {
  AbstractRefinement::refinement_list_type list;

  // this refinement algorithm uses the SVM refinement indicator.
  // dynamic casting is used to maintain the signature of the algorithm,
  // but still be able to use the
  // SVM refinement indicator with it.
  const ForwardSelectorRefinementIndicator& svmIndicator =
    dynamic_cast<const ForwardSelectorRefinementIndicator&>(functor);
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
    
    //std::cout << "refinable point : " << storage.getSequenceNumber(point) << std::endl;
    double measure = 1.0 / svmIndicator(storage, storage.getSequenceNumber(point));

    if (measure > iThreshold_) {
      size_t d = 0;

      key = new refinement_key_type(point, storage.getSequenceNumber(point), d);        // d actually not required now
      list.emplace_front(std::shared_ptr<AbstractRefinement::refinement_key_type>(key),
                         measure); 
    }
  }

  return list;
}


void ForwardSelectorRefinement::addElementToCollection(
  const GridStorage::grid_map_iterator& iter,
  AbstractRefinement::refinement_list_type current_value_list,
  size_t refinements_num,
  AbstractRefinement::refinement_container_type& collection) {
  for (AbstractRefinement::refinement_list_type::iterator it =
         current_value_list.begin(); it != current_value_list.end(); it++) {
    collection.push_back(*it);
    std::push_heap(collection.begin(), collection.end(), AbstractRefinement::compare_pairs);

    if (collection.size() > refinements_num) {
      // remove the top (biggest) element
      std::pop_heap(collection.begin(), collection.end(),
                    AbstractRefinement::compare_pairs);
      collection.pop_back();
    }

  }
}


void ForwardSelectorRefinement::refineGridpointsCollection(
  GridStorage& storage, RefinementFunctor& functor,
  AbstractRefinement::refinement_container_type& collection) {
  //ForwardSelectorRefinementIndicator::value_type max_value;

  // now refine all grid points which satisfy the refinement criteria

  //double threshold = functor.getRefinementThreshold();
  refinement_key_type* key;

  ForwardSelectorRefinementIndicator& svmIndicator =
    dynamic_cast<ForwardSelectorRefinementIndicator&>(functor);

  for (AbstractRefinement::refinement_pair_type& pair : collection) {
    key = dynamic_cast<refinement_key_type*>(pair.first.get());

    // initialize new components of w1 and w2 using insertable grid points
    //ForwardSelectorRefinement::insertables_container_type insertables; // better: shared pointer to insertables-container
    GridPoint& point = key->getPoint();
    /*GridStorage::grid_map_iterator child_iter;
    GridStorage::grid_map_iterator end_iter = storage.end();*/

    // collect insertable grid points (new children) in dimension key->getDim() -> no for loop
    // maybe better to collect all insertables for a point at once and check if a grid point has already been seen
     
    //const ForwardSelectorRefinementIndicator& svmIndicator =
    //  dynamic_cast<const ForwardSelectorRefinementIndicator&>(functor);
    /*ForwardSelectorRefinementIndicator& svmIndicator =
      dynamic_cast<ForwardSelectorRefinementIndicator&>(functor);*/

    //for (size_t d = 0; d < storage.getDimension(); d++) {
    /*size_t d = key->getDim();
    index_t source_index;
    level_t source_level;
    point.get(d, source_level, source_index);

    // test existence of left child
    point.set(d, source_level + 1, 2 * source_index - 1);
    child_iter = storage.find(&point);

    if (child_iter == end_iter) {
      //insertables.push_back(point);
      svmIndicator.update(point);
    }

    // test existance of right child
    point.set(d, source_level + 1, 2 * source_index + 1);
    child_iter = storage.find(&point);

    if (child_iter == end_iter) {
      //insertables.push_back(point); // point or &point ?
      svmIndicator.update(point);
    }

    // reset current grid point in dimension d
    point.set(d, source_level, source_index);*/

    //}
    // compute new w1 and w2 components
    
    //svmIndicator.update(insertables);

    // check if insertables benefit goal function

    // refine grid point
    //if (pair.second > functor.start() && pair.second >= threshold) {
    //}

    //std::cout << "refine point: " << storage.getSequenceNumber(point) << std::endl;

    //this->refineGridpoint1D(storage, key->getPoint(), key->getDim());
    //this->refineGridpoint(storage, storage.getSequenceNumber(point), );
    //refineGridpointSVM(storage, storage.getSequenceNumber(point), svmIndicator);

    storage[storage.getSequenceNumber(point)].setLeaf(false);
    //point.setLeaf(false);
    
    //void ForwardSelectorRefinementIndicator::*update_w = &ForwardSelectorRefinementIndicator::update;

    for (size_t d = 0; d < storage.getDimension(); d++) {
      index_t source_index;
      level_t source_level;
      point.get(d, source_level, source_index);
      // generate left child, if necessary
      point.set(d, source_level + 1, 2 * source_index - 1);

      if (!storage.isContaining(point)) {
        point.setLeaf(true);
        svmCreateGridpoint(storage, point, svmIndicator);
      }

      // generate right child, if necessary
      point.set(d, source_level + 1, 2 * source_index + 1);
 
      if (!storage.isContaining(point)) {
        point.setLeaf(true);
        svmCreateGridpoint(storage, point, svmIndicator);
      }
 
      point.set(d, source_level, source_index);
    }
    
  }

  collection.empty();
}


void ForwardSelectorRefinement::svmCreateGridpoint(
  GridStorage& storage, GridPoint& point, ForwardSelectorRefinementIndicator& svmIndicator) {
  index_t source_index;
  level_t source_level;

  for (size_t d = 0; d < storage.getDimension(); d++) {
    svmCreateGridpoint1D(point, d, storage, source_index, source_level, svmIndicator);
  }

  storage.insert(point);
  svmIndicator.update(point);
}


void ForwardSelectorRefinement::svmCreateGridpoint1D(GridPoint& point, size_t d, GridStorage& storage,
                                                     index_t& source_index, level_t& source_level,
                                                     ForwardSelectorRefinementIndicator& svmIndicator) {
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
      svmCreateGridpoint(storage, point, svmIndicator);
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


void ForwardSelectorRefinement::free_refine(GridStorage& storage, ForwardSelectorRefinementIndicator& functor) {
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
