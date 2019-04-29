// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/generation/hashmap/HashCoarsening.hpp>

#include <sgpp/globaldef.hpp>

#include <cmath>
#include <list>
#include <utility>
#include <vector>
#include <algorithm>
#include <chrono>
#include <ctime> 
namespace sgpp {
namespace base {

void HashCoarsening::free_coarsen_NFirstOnly(GridStorage& storage,
                                             CoarseningFunctor& functor,
                                             DataVector& alpha,
                                             size_t numFirstPoints,
                                             size_t minIndexConsidered,
                                             std::vector<HashGridPoint>* removedPoints,
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
    std::cout<<"(";
    for (size_t d = 0; d < point.getDimension(); d++){
        std::cout<<point.getStandardCoordinate(d);
        if(d!=point.getDimension()-1){
            std::cout<<",";
        }
        else{
            std::cout<<")";
        }
    }
    if (point.isLeaf()) {
      auto start = std::chrono::system_clock::now();

      CoarseningFunctor::value_type current_value = functor(storage, z);
          
      auto end = std::chrono::system_clock::now();
      std::chrono::duration<double> elapsed_seconds = end-start;
      std::cout<< "++++++++++++++++elapsed time: " << elapsed_seconds.count() << "s\n";
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
    else{

        std::cout<<";NaN;NaN;NaN"<<std::endl;
    }
  }

  // DEBUG : print list of removable candidates
  std::cout << "list of removable candidates:\n";
  for (size_t i = 0; i < remove_num; i++) {
    std::cout << "Index: " << removeCandidates[i].first << " with surplus " <<
    removeCandidates[i].second
              << std::endl;
  }
  std::cout << std::endl;

  // remove the marked grid point if their surplus
  // is below the given threshold
  CoarseningFunctor::value_type threshold = functor.getCoarseningThreshold();
  CoarseningFunctor::value_type initValue = functor.start();

  // vector to save remaining points
  std::vector<size_t> remainingIndex;

  // vector containing the actually removed points. "local" b/c @param removedPoints
  std::vector<size_t> localRemovedPoints;

  // vector to store the points that match all condition for deleting
  // this->deletePoints.clear();

  for (size_t i = 0; i < remove_num; i++) {
    if (removeCandidates[i].second < initValue && removeCandidates[i].second <= threshold) {
      localRemovedPoints.push_back(removeCandidates[i].first);
      std::cout<<"This candidate has score:"<<removeCandidates[i].second<<std::endl;
      if (removedPoints != 0) {
        removedPoints->push_back(GridPoint(storage.getPoint(removeCandidates[i].first)));
      }
      if (removedSeq != 0) {
        removedSeq->push_back(removeCandidates[i].first);
      }
    }
  }

  // DEBUG : print list points to delete



  // For some reason HashGridStorage expects a std::list and not a vector D:
  std::list<size_t> removedPointsList(localRemovedPoints.begin(),
                                      localRemovedPoints.end());
  std::cout << "list of points to delete:\n";
  for (auto i:removedPointsList){
      std::cout << i << ",";
  }
  std::cout<<std::endl;
  remainingIndex = storage.deletePoints(removedPointsList);

  // DEBUG
  // std::cout << "List of remaining GridPoints (indices)" << std::endl;
  // for (size_t i = 0; i < remainingIndex.size(); i++)
  // {
  //   std::cout << remainingIndex[i] << " ";
  // }
  // std::cout << std::endl << std::endl;

  // Drop Elements from DataVector
  alpha.restructure(remainingIndex);

  delete[] removeCandidates;
}

void HashCoarsening::free_coarsen(GridStorage& storage,
                                  CoarseningFunctor& functor,
                                  DataVector& alpha,
                                  std::vector<HashGridPoint>* removedPoints,
                                  std::vector<size_t>* removedSeq) {
  free_coarsen_NFirstOnly(storage, functor, alpha, storage.getSize(), 0, removedPoints, removedSeq);
}


size_t HashCoarsening::getNumberOfRemovablePoints(GridStorage& storage) {
  size_t counter = 0;

  if (storage.getSize() == 0) {
    throw generation_exception("storage empty");
  }

  GridPoint point;
  GridStorage::grid_map_iterator end_iter = storage.end();

  for (GridStorage::grid_map_iterator iter = storage.begin(); iter != end_iter; iter++) {
    point = *(iter->first);

    if (point.isLeaf()) {
      counter++;
    }
  }

  return counter;
}

}  // namespace base
}  // namespace sgpp
