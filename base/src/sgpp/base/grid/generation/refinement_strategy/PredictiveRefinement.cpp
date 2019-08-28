// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/generation/refinement_strategy/PredictiveRefinement.hpp>

#include <sgpp/globaldef.hpp>
#include <sgpp/base/grid/generation/functors/PredictiveRefinementIndicator.hpp>

#include <algorithm>
#include <limits>
#include <list>
#include <vector>
#include <utility>
#include <iterator>


namespace sgpp {
namespace base {

bool doubleReverseCompare(const double firstEl, const double secondEl) {
  return firstEl > secondEl;
}


void PredictiveRefinement::addElementToCollection(
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


AbstractRefinement::refinement_list_type PredictiveRefinement::getIndicator(
  GridStorage& storage,
  const GridStorage::grid_map_iterator& iter,
  const RefinementFunctor& functor) const {
  AbstractRefinement::refinement_list_type list;

  // this refinement algorithm uses the predictive refinement indicator.
  // dynamic casting is used to maintain the signature of the algorithm,
  // but still be able to use the
  // predictive refinement indicator with it.
  const PredictiveRefinementIndicator& errorIndicator =
    dynamic_cast<const PredictiveRefinementIndicator&>(functor);
  refinement_key_type* key;

  GridPoint& point = *(iter->first);
  GridStorage::grid_map_iterator child_iter;
  GridStorage::grid_map_iterator end_iter = storage.end();

  for (size_t d = 0; d < storage.getDimension(); d++) {
    index_t source_index;
    level_t source_level;
    point.get(d, source_level, source_index);
    double error = errorIndicator.start();
    // errorIndicator->setActiveDim(d);

    // test existence of left child
    point.set(d, source_level + 1, 2 * source_index - 1);
    child_iter = storage.find(&point);

    if (child_iter == end_iter) {
      // use the predictive error indicator
      error += errorIndicator(point);
    }

    // test existance of right child
    point.set(d, source_level + 1, 2 * source_index + 1);
    child_iter = storage.find(&point);

    if (child_iter == end_iter) {
      // use predictive refinement indicator
      // use the predictive error indicator
      error += errorIndicator(point);
    }

    // reset current grid point in dimension d
    point.set(d, source_level, source_index);

    if (error > iThreshold_) {
      key = new refinement_key_type(*(iter->first),
                                    storage.getSequenceNumber(*iter->first), d);
      list.emplace_front(
        std::shared_ptr<AbstractRefinement::refinement_key_type>(key),
        error);
    }
  }

  return list;
}


void PredictiveRefinement::collectRefinablePoints(
  GridStorage& storage, RefinementFunctor& functor,
  AbstractRefinement::refinement_container_type& collection) {
  size_t refinements_num = functor.getRefinementsNum();

  GridPoint index;
  GridStorage::grid_map_iterator end_iter = storage.end();

  // start iterating over whole grid
  for (GridStorage::grid_map_iterator iter = storage.begin();
       iter != end_iter; iter++) {
    index = *(iter->first);
    AbstractRefinement::refinement_list_type current_value_list = getIndicator(
          storage, iter, functor);
    addElementToCollection(iter, current_value_list, refinements_num,
                           collection);
  }
}


void PredictiveRefinement::refineGridpointsCollection(
  GridStorage& storage, RefinementFunctor& functor,
  AbstractRefinement::refinement_container_type& collection) {
  // now refine all grid points which satisfy the refinement criteria
  double threshold = functor.getRefinementThreshold();
  refinement_key_type* key;

  for (AbstractRefinement::refinement_pair_type& pair : collection) {
    key = dynamic_cast<refinement_key_type*>(pair.first.get());

    if (pair.second > functor.start() && pair.second >= threshold) {
      this->refineGridpoint1D(storage, key->getPoint(), key->getDim());
      key->getPoint().setLeaf(false);
    }

    // delete key;
  }

  collection.empty();
}

void PredictiveRefinement::free_refine(GridStorage& storage,
                                       PredictiveRefinementIndicator& functor,
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

  // initialization
  AbstractRefinement::refinement_container_type collection;

  collectRefinablePoints(storage, functor, collection);
  // now refine all grid points which satisfy the refinement criteria
  refineGridpointsCollection(storage, functor, collection);
  collection.clear();

  if (addedPoints != nullptr) {
    for (size_t i = sizeBeforeRefine; i < storage.getSize(); i++) {
      addedPoints->push_back(i);
    }
  }
}

}  // namespace base
}  // namespace sgpp
