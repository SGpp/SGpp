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

  double threshold = svmIndicator.getRefinementThreshold();
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

    // if child does not exist -> check indicator
    if (child_iter == storage.end()) {
      isRefinablePoint = true;
      point.set(d, source_level, source_index);
      break;
    }

    // test existance of right child
    point.set(d, source_level + 1, 2 * source_index + 1);
    child_iter = storage.find(&point);

    // if child does not exist -> check indicator
    if (child_iter == storage.end()) {
      isRefinablePoint = true;
      point.set(d, source_level, source_index);
      break;
    }

    // reset current grid point in dimension d
    point.set(d, source_level, source_index);
  }

  if (isRefinablePoint) {
    
    //std::cout << "refinable point: " << storage.getSequenceNumber(point) << std::endl;
    double measure = 1.0 / svmIndicator(storage, storage.getSequenceNumber(point));

    if (measure > threshold) {
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

  // now refine all grid points which satisfy the refinement criteria
  refinement_key_type* key;

  ForwardSelectorRefinementIndicator& svmIndicator =
    dynamic_cast<ForwardSelectorRefinementIndicator&>(functor);

  // check last sequence number
  size_t lastSeqNr = storage.getSize() - 1;

  for (AbstractRefinement::refinement_pair_type& pair : collection) {
    key = dynamic_cast<refinement_key_type*>(pair.first.get());
    GridPoint& point = key->getPoint();
    this->refineGridpoint(storage, storage.getSequenceNumber(point));
    //point.setLeaf(false); // this is done within refineGridpoint() already    
  }
  // extend w1 and w2 vectors
  for (size_t seqNr = lastSeqNr + 1; seqNr < storage.getSize(); ++seqNr) {
    svmIndicator.update(storage.getPoint(seqNr));
  }
  collection.empty();
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
