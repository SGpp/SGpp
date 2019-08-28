// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/generation/hashmap/HashRefinementBoundaries.hpp>
#include <sgpp/base/exception/generation_exception.hpp>

#include <sgpp/globaldef.hpp>

#include <cmath>
#include <vector>
#include <memory>


namespace sgpp {
namespace base {


void HashRefinementBoundaries::addElementToCollection(
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


AbstractRefinement::refinement_list_type HashRefinementBoundaries::getIndicator(
  GridStorage& storage,
  const GridStorage::grid_map_iterator& iter,
  const RefinementFunctor& functor) const {
  AbstractRefinement::refinement_list_type list;
  list.emplace_front(std::make_shared<AbstractRefinement::refinement_key_type>(*
                     (iter->first), iter->second),
                     functor(storage, iter->second));
  return list;
}



void HashRefinementBoundaries::collectRefinablePoints(GridStorage& storage,
    RefinementFunctor& functor,
    AbstractRefinement::refinement_container_type& collection) {

  size_t refinements_num = functor.getRefinementsNum();
  GridPoint point;
  GridStorage::grid_map_iterator end_iter = storage.end();

  // I think this may be dependent on local support
  for (GridStorage::grid_map_iterator iter = storage.begin(); iter != end_iter;
       iter++) {
    point = *(iter->first);

    GridStorage::grid_map_iterator child_iter;

    for (size_t d = 0; d < storage.getDimension(); d++) {
      index_t source_index;
      level_t source_level;
      point.get(d, source_level, source_index);

      if (source_level == 0) {
        // we only have one child on level 1
        point.set(d, 1, 1);
        child_iter = storage.find(&point);

        // if there no more grid points --> test if we should refine the grid
        if (child_iter == end_iter) {
          AbstractRefinement::refinement_list_type current_value_list =
            getIndicator(storage, iter, functor);
          addElementToCollection(iter, current_value_list, refinements_num,
                                 collection);
          break;
        }
      } else {
        // left child
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

        // right child
        point.set(d, source_level + 1, 2 * source_index + 1);
        child_iter = storage.find(&point);

        if (child_iter == end_iter) {
          AbstractRefinement::refinement_list_type current_value_list =
            getIndicator(storage, iter, functor);
          addElementToCollection(iter, current_value_list, refinements_num,
                                 collection);
          break;
        }
      }

      point.set(d, source_level, source_index);
    }
  }
}


void HashRefinementBoundaries::refineGridpointsCollection(GridStorage& storage,
    RefinementFunctor& functor,
    AbstractRefinement::refinement_container_type& collection) {
  double threshold = functor.getRefinementThreshold();

  for (AbstractRefinement::refinement_pair_type& pair : collection) {
    if (pair.second >= threshold) {
      refineGridpoint(storage, pair.first->getSeq());
    }
  }
}

void HashRefinementBoundaries::free_refine(GridStorage& storage,
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
  // can refine grid on several points
  refineGridpointsCollection(storage, functor, collection);

  if (addedPoints != nullptr) {
    for (size_t i = sizeBeforeRefine; i < storage.getSize(); i++) {
      addedPoints->push_back(i);
    }
  }
}


size_t HashRefinementBoundaries::getNumberOfRefinablePoints(
  GridStorage& storage) {
  size_t counter = 0;

  if (storage.getSize() == 0) {
    throw generation_exception("storage empty");
  }

  GridPoint point;
  GridStorage::grid_map_iterator end_iter = storage.end();

  // I think this may be dependent on local support
  for (GridStorage::grid_map_iterator iter = storage.begin(); iter != end_iter;
       iter++) {
    point = *(iter->first);

    GridStorage::grid_map_iterator child_iter;

    for (size_t d = 0; d < storage.getDimension(); d++) {
      index_t source_index;
      level_t source_level;
      point.get(d, source_level, source_index);

      if (source_level == 0) {
        // level 1
        point.set(d, 1, 1);
        child_iter = storage.find(&point);

        // if there no more grid points --> test if we should refine the grid
        if (child_iter == end_iter) {
          counter++;
          break;
        }
      } else {
        // left child
        point.set(d, source_level + 1, 2 * source_index - 1);
        child_iter = storage.find(&point);

        // if there no more grid points --> test if we should refine the grid
        if (child_iter == end_iter) {
          counter++;
          break;
        }

        // right child
        point.set(d, source_level + 1, 2 * source_index + 1);
        child_iter = storage.find(&point);

        if (child_iter == end_iter) {
          counter++;
          break;
        }
      }

      point.set(d, source_level, source_index);
    }
  }

  return counter;
}


void HashRefinementBoundaries::refineGridpoint1D(GridStorage& storage,
    GridPoint& point, size_t d) {
  index_t source_index;
  level_t source_level;
  point.get(d, source_level, source_index);

  if (source_level == 0) {
    // we only have one child on level 1
    point.set(d, 1, 1);

    if (!storage.isContaining(point)) {
      point.setLeaf(true);
      createGridpoint(storage, point);
    }
  } else {
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
  }

  point.set(d, source_level, source_index);
}


void HashRefinementBoundaries::refineGridpoint(GridStorage& storage,
    size_t refine_index) {
  GridPoint point(storage[refine_index]);

  // Sets leaf property of index, which is refined to false
  storage[refine_index].setLeaf(false);

  for (size_t d = 0; d < storage.getDimension(); d++) {
    refineGridpoint1D(storage, point, d);
  }
}



void HashRefinementBoundaries::createGridpoint(GridStorage& storage,
    GridPoint& point) {
  // create grid with its needed childern and parents
  createGridpointGeneral(storage, point);
  // create all missing points an level zero
  createGridpointLevelZeroConsistency(storage, point);
}

void HashRefinementBoundaries::createGridpoint1D(GridPoint& point,
    size_t d, GridStorage& storage, index_t& source_index,
    level_t& source_level) {
  point.get(d, source_level, source_index);

  if (source_level == 1) {
    // check if we need some additional points on the boundaries,
    // only needed on a N dim grid
    if (storage.getDimension() > 1) {
      // test if there are boundaries in every dimension for this grid point
      // left boundary
      point.set(d, 0, 0);
      createGridpointSubroutine(storage, point);

      // right boundary
      point.set(d, 0, 1);
      createGridpointSubroutine(storage, point);

      // restore values
      point.set(d, source_level, source_index);
    }
  }

  AbstractRefinement::createGridpoint1D(point, d, storage, source_index,
                                        source_level);
}

void HashRefinementBoundaries::createGridpointGeneral(GridStorage& storage,
    GridPoint& point) {
  index_t source_index;
  level_t source_level;

  for (size_t d = 0; d < storage.getDimension(); d++) {
    createGridpoint1D(point, d, storage, source_index, source_level);
  }

  storage.insert(point);
}


void HashRefinementBoundaries::createGridpointLevelZeroConsistency(
  GridStorage& storage, GridPoint& point) {
  for (size_t d = 0; d < storage.getDimension(); d++) {
    index_t source_index;
    level_t source_level;
    point.get(d, source_level, source_index);

    // Assure that we have always a consistent grid with both functions
    // 0,0 and 0,1 on level zero
    if (source_level == 0) {
      // check if we need some additional points on the boundaries,
      // only needed on a N dim grid
      if (storage.getDimension() > 1) {
        // if we have already a left boundary...
        point.set(d, 0, 0);

        if (storage.isContaining(point)) {
          // ... we have to read leaf property
          bool Leaf = point.isLeaf();
          // ... we have to generate the correspondending right boundary
          point.set(d, 0, 1);

          if (!storage.isContaining(point)) {
            bool saveLeaf = point.isLeaf();
            point.setLeaf(Leaf);
            createGridpoint(storage, point);
            point.setLeaf(saveLeaf);
          } else {
            // set stored index to Leaf from the left boundary
            storage.getPoint(storage.find(&point)->second).setLeaf(Leaf);
          }
        }

        // if we have already a right boundary...
        point.set(d, 0, 1);

        if (storage.isContaining(point)) {
          // ... we have to read leaf property
          bool Leaf = point.isLeaf();
          // ... we have to generate the correspondending right boundary
          point.set(d, 0, 0);

          if (!storage.isContaining(point)) {
            bool saveLeaf = point.isLeaf();
            point.setLeaf(Leaf);
            createGridpoint(storage, point);
            point.setLeaf(saveLeaf);
          } else {
            // set stored index to Leaf from the right boundary
            storage.getPoint(storage.find(&point)->second).setLeaf(Leaf);
          }
        }

        // restore values
        point.set(d, source_level, source_index);
      }
    }
  }
}

}  // namespace base
}  // namespace sgpp
