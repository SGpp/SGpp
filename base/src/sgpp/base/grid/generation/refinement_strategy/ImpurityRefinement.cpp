// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/base/grid/generation/refinement_strategy/ImpurityRefinement.hpp>
#include <sgpp/base/grid/generation/functors/ImpurityRefinementIndicator.hpp>

#include <list>
#include <vector>
#include <iterator>


namespace sgpp {
namespace base {


void ImpurityRefinement::collectRefinablePoints(GridStorage& storage, RefinementFunctor& functor, 
  AbstractRefinement::refinement_container_type& collection) {

  size_t refinementsNum = functor.getRefinementsNum();
  GridStorage::grid_map_iterator end_iter = storage.end();

  // start iterating over whole grid
  for (GridStorage::grid_map_iterator iter = storage.begin(); iter != end_iter; iter++) {
    AbstractRefinement::refinement_list_type current_value_list = 
      getIndicator(storage, iter, functor);
    addElementToCollection(iter, current_value_list, refinementsNum, collection);
  }
}

AbstractRefinement::refinement_list_type ImpurityRefinement::getIndicator(
  GridStorage& storage,
  const GridStorage::grid_map_iterator& iter,
  const RefinementFunctor& functor) const {
  AbstractRefinement::refinement_list_type list;

  // this refinement algorithm uses the impurity refinement indicator
  const ImpurityRefinementIndicator& impurityIndicator =
    dynamic_cast<const ImpurityRefinementIndicator&>(functor);
  refinement_key_type* key;

  GridPoint& point = *(iter->first);
  GridStorage::grid_map_iterator child_iter;
  GridStorage::grid_map_iterator end_iter = storage.end();

  double threshold = impurityIndicator.getRefinementThreshold();
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

    // if child is missing,
    // test if we should refine the grid
    if (child_iter == storage.end()) {
      isRefinablePoint = true;
      point.set(d, source_level, source_index);
      break;
    }

    // test existance of right child
    point.set(d, source_level + 1, 2 * source_index + 1);
    child_iter = storage.find(&point);

    // if child is missing,
    // test if we should refine the grid
    if (child_iter == storage.end()) {
      isRefinablePoint = true;
      point.set(d, source_level, source_index);
      break;
    }

    // reset current grid point in dimension d
    point.set(d, source_level, source_index);
  }

  if (isRefinablePoint) {
    // evaluate indicator
    double impurity = impurityIndicator(point);

    if (impurity > threshold) {
      size_t d = 0;
      key = new refinement_key_type(point, storage.getSequenceNumber(point), d);        
      list.emplace_front(std::shared_ptr<AbstractRefinement::refinement_key_type>(key),
                         impurity); 
    }
  }

  return list;
}


void ImpurityRefinement::addElementToCollection(
  const GridStorage::grid_map_iterator& iter,
  AbstractRefinement::refinement_list_type current_value_list,
  size_t refinementsNum,
  AbstractRefinement::refinement_container_type& collection) {

  for (AbstractRefinement::refinement_list_type::iterator it =
         current_value_list.begin(); it != current_value_list.end(); it++) {
    collection.push_back(*it);
    std::push_heap(collection.begin(), collection.end(), AbstractRefinement::compare_pairs);

    if (collection.size() > refinementsNum) {
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

  refinement_key_type* key;

  ImpurityRefinementIndicator& impurityIndicator =
    dynamic_cast<ImpurityRefinementIndicator&>(functor);

  // check last sequence number
  size_t lastSeqNr = storage.getSize() - 1;

  for (AbstractRefinement::refinement_pair_type& pair : collection) {
    key = dynamic_cast<refinement_key_type*>(pair.first.get());

    GridPoint& point = key->getPoint();

    //std::cout << "refine point: " << storage.getSequenceNumber(point) << std::endl;
    this->refineGridpoint(storage, storage.getSequenceNumber(point));

    point.setLeaf(false);     
  }
  // for SVM learner -> extend w1 and w2 vectors
  if ( (impurityIndicator.w1 != nullptr) 
  &&   (impurityIndicator.w1 != nullptr) 
  &&   (impurityIndicator.alphas != nullptr) ) {
    for (size_t seqNr = lastSeqNr + 1; seqNr < storage.getSize(); ++seqNr) {
      impurityIndicator.update(storage.getPoint(seqNr));
    }
  }
  collection.empty();
}


void ImpurityRefinement::free_refine(GridStorage& storage, 
                                     ImpurityRefinementIndicator& functor) {
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
