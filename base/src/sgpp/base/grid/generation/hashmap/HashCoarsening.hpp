// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef HASHCOARSENING_HPP
#define HASHCOARSENING_HPP

#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/generation/functors/CoarseningFunctor.hpp>

#include <sgpp/base/exception/generation_exception.hpp>

#include <sgpp/globaldef.hpp>

#include <vector>
#include <list>
#include <cmath>
#include <utility>


namespace SGPP {
namespace base {

/**
 * Standard free coarsening class for sparse grids, only
 * inner grid points can be removed
 */
class HashCoarsening {
 public:
  typedef GridStorage::index_type index_type;
  typedef index_type::index_type index_t;
  typedef index_type::level_type level_t;
  typedef std::pair<size_t, CoarseningFunctor::value_type> GridPoint;

  /**
   * Performs coarsening on grid. It's possible to remove a certain number
   * of gridpoints in one coarsening step. This number is specified within the
   * declaration of the coarsening functor. Also the coarsening threshold is
   * specified in the coarsening functor. ONLY INNER GRID POINTS WILL
   * BE REMOVED!
   *
   * Here only the numFirstPoints are regarded for coarsening, later points
   * are skipped.
   *
   * @param storage hashmap that stores the grid points
   * @param functor a function used to determine if refinement is needed
   * @param alpha pointer to the gridpoints' coefficients removed points must also be considered in this vector
   * @param numFirstPoints number of grid points that are regarded to be coarsened
   */
  void free_coarsen_NFirstOnly(GridStorage& storage, CoarseningFunctor& functor,
                               DataVector& alpha, size_t numFirstPoints) {
    // check if the grid has any points
    if (storage.size() == 0) {
      throw generation_exception("storage empty");
    }

    // Perepare temp-data in order to determine the removable grid points
    // -> leafs with minimal surplus
    size_t remove_num = functor.getRemovementsNum();

    if (remove_num == 0) return;

    // create an array that will contain the GridPoints
    // (pair of the grid Point's index and its surplus)
    // that should be removed
    GridPoint* removePoints = new GridPoint[remove_num];

    // init the removePoints array:
    // set initial surplus and set all indices to zero
    for (size_t i = 0; i < remove_num; i++) {
      removePoints[i].second = functor.start();
      removePoints[i].first = 0;
    }

    // help variable to store the gridpoint with highest
    // surplus in removePoints
    size_t max_idx = 0;

    // assure that only the first numFirstPoints are checked for coarsening
    for (size_t z = 0; z < numFirstPoints; z++) {
      index_type* index = storage.get(z);

      if (index->isLeaf() && index->isInnerPoint()) {
        CoarseningFunctor::value_type current_value = functor(storage, z);

        if (current_value < removePoints[max_idx].second) {
          // Replace the maximum point array of removable candidates,
          // find the new maximal point
          removePoints[max_idx].second = current_value;
          removePoints[max_idx].first = z;

          // find new maximum entry
          max_idx = 0;

          for (size_t i = 1; i < remove_num; i++) {
            if (removePoints[i].second > removePoints[max_idx].second) {
              max_idx = i;
            }
          }
        }
      }
    }

    // DEBUG : print list of removable candidates
    // for (size_t i = 0; i < remove_num; i++)
    // {
    //   std::cout << "Index: " << removePoints[i].first <<
    //   " with surplus " << removePoints[i].second << std::endl;
    // }
    // std::cout << std::endl;

    // remove the marked grid point if their surplus
    // is below the given threshold
    CoarseningFunctor::value_type threshold = functor.getCoarseningThreshold();
    CoarseningFunctor::value_type initValue = functor.start();

    // vector to save remaining points
    std::vector<size_t> remainingIndex;

    // vector to stored the points that match all condition for deleting
    std::list<size_t> deletePoints;

    for (size_t i = 0; i < remove_num; i++) {
      if (removePoints[i].second < initValue &&
          removePoints[i].second <= threshold) {
        deletePoints.push_back(removePoints[i].first);
      }
    }

    // DEBUG : print list points to delete
    // for(std::list<size_t>::iterator iter = deletePoints.begin();
    // iter != deletePoints.end(); iter++)
    // {
    //   std::cout << "Index: " << *iter << std::endl;
    // }

    remainingIndex = storage.deletePoints(deletePoints);

    // DEBUG
    // std::cout << "List of remaining GridPoints (indices)" << std::endl;
    // for (size_t i = 0; i < remainingIndex.size(); i++)
    // {
    //   std::cout << remainingIndex[i] << " ";
    // }
    // std::cout << std::endl << std::endl;

    // Drop Elements from DataVector
    alpha.restructure(remainingIndex);

    delete[] removePoints;
  }

  /**
   * Performs coarsening on grid. It's possible to remove a certain number
   * of gridpoints in one coarsening step. This number is specified within the
   * declaration of the coarsening functor. Also the coarsening threshold is
   * specified in the coarsening functor. ONLY INNER GRID POINTS WILL
   * BE REMOVED!
   *
   * This function calls free_coarsen_NFirstOnly with numFirstPoints equal
   * to the grid's size.
   *
   * @param storage hashmap that stores the grid points
   * @param functor a function used to determine if refinement is needed
   * @param alpha pointer to the gridpoints' coefficients removed points must also be considered in this vector
   */
  void free_coarsen(GridStorage& storage, CoarseningFunctor& functor,
                    DataVector& alpha) {
    free_coarsen_NFirstOnly(storage, functor, alpha, storage.size());
  }

  /**
   * Calculates the number of points, which can be refined
   *
   * @param storage hashmap that stores the grid points
   */
  size_t getNumberOfRemovablePoints(GridStorage& storage) {
    size_t counter = 0;

    if (storage.size() == 0) {
      throw generation_exception("storage empty");
    }

    index_type index;
    GridStorage::grid_map_iterator end_iter = storage.end();

    for (GridStorage::grid_map_iterator iter = storage.begin();
         iter != end_iter;
         iter++) {
      index = *(iter->first);

      if (index.isLeaf()) {
        counter++;
      }
    }

    return counter;
  }
};

}  // namespace base
}  // namespace SGPP

#endif /* HASHCOARSENING_HPP */
