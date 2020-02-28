// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/generation/hashmap/HashCoarsening.hpp>
#include <sgpp/globaldef.hpp>

#include <algorithm>
#include <cmath>
#include <list>
#include <utility>
#include <vector>

namespace sgpp {
namespace base {

void HashCoarsening::free_coarsen_NFirstOnly(
    GridStorage& storage, CoarseningFunctor& functor, size_t numFirstPoints,
    size_t minIndexConsidered, std::vector<HashGridPoint>* removedPoints,
    std::vector<size_t>* removedSeq) {
  // check if the grid has any points
  if (storage.getSize() == 0) {
    throw generation_exception("storage empty");
  }

  // Perepare temp-data in order to determine the removable grid points
  // -> leafs with minimal surplus
  // Makes sure at most grid-size number of points are considered for coarsening
  size_t remove_num = std::min(storage.getSize(), functor.getRemovementsNum());

  if (remove_num == 0) return;

  // create an array that will contain the GridPoints
  // (pair of the grid Point's index and its surplus)
  // that should be removed
  typedef std::pair<size_t, CoarseningFunctor::value_type> GridPointPair;
  GridPointPair* removeCandidates = new GridPointPair[remove_num];

  // init the removeCandidates array:
  // set initial surplus and set all indices to zero
  for (size_t i = 0; i < remove_num; i++) {
    removeCandidates[i].second = functor.start();
    removeCandidates[i].first = 0;
  }

  // help variable to store the gridpoint with highest
  // surplus in removeCandidates
  size_t max_idx = 0;

  // assure that only the first numFirstPoints are checked for coarsening
  // also assure, that indices bigger than minIndexConsidered are not checked
  for (size_t z = minIndexConsidered; z < numFirstPoints; z++) {
    GridPoint& point = storage.getPoint(z);

    if (point.isLeaf()) {
      CoarseningFunctor::value_type current_value = functor(storage, z);

      if (current_value < removeCandidates[max_idx].second) {
        // Replace the maximum point array of removable candidates,
        // find the new maximal point
        removeCandidates[max_idx].second = current_value;
        removeCandidates[max_idx].first = z;

        // find new maximum entry
        max_idx = 0;

        for (size_t i = 1; i < remove_num; i++) {
          if (removeCandidates[i].second > removeCandidates[max_idx].second) {
            max_idx = i;
          }
        }
      }
    }
  }

  // remove the marked grid point if their surplus
  // is below the given threshold
  CoarseningFunctor::value_type threshold = functor.getCoarseningThreshold();
  CoarseningFunctor::value_type initValue = functor.start();

  // vector to save remaining points
  std::vector<size_t> remainingIndex;

  // vector containing the actually removed points. "local" b/c @param
  // removedPoints
  std::vector<size_t> localRemovedPoints;

  // vector to store the points that match all condition for deleting
  // this->deletePoints.clear();

  for (size_t i = 0; i < remove_num; i++) {
    if (removeCandidates[i].second < initValue &&
        removeCandidates[i].second <= threshold) {
      localRemovedPoints.push_back(removeCandidates[i].first);
      if (removedPoints != nullptr) {
        removedPoints->push_back(
            GridPoint(storage.getPoint(removeCandidates[i].first)));
      }
      if (removedSeq != nullptr) {
        removedSeq->push_back(removeCandidates[i].first);
      }
    }
  }

  // For some reason HashGridStorage expects a std::list and not a vector D:
  std::list<size_t> removedPointsList(localRemovedPoints.begin(),
                                      localRemovedPoints.end());
  storage.deletePoints(removedPointsList);

  delete[] removeCandidates;
}

void HashCoarsening::free_coarsen(GridStorage& storage,
                                  CoarseningFunctor& functor,
                                  std::vector<HashGridPoint>* removedPoints,
                                  std::vector<size_t>* removedSeq) {
  free_coarsen_NFirstOnly(storage, functor, storage.getSize(), 0, removedPoints,
                          removedSeq);
}

size_t HashCoarsening::getNumberOfRemovablePoints(GridStorage& storage) {
  size_t counter = 0;

  if (storage.getSize() == 0) {
    throw generation_exception("storage empty");
  }

  GridPoint point;
  GridStorage::grid_map_iterator end_iter = storage.end();

  for (GridStorage::grid_map_iterator iter = storage.begin(); iter != end_iter;
       iter++) {
    point = *(iter->first);

    if (point.isLeaf()) {
      counter++;
    }
  }

  return counter;
}

}  // namespace base
}  // namespace sgpp
