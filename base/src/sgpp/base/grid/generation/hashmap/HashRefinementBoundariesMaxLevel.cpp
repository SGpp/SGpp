// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/generation/hashmap/HashRefinementBoundariesMaxLevel.hpp>
#include <sgpp/base/exception/generation_exception.hpp>

#include <sgpp/globaldef.hpp>

#include <vector>
#include <cmath>


namespace sgpp {
namespace base {

void HashRefinementBoundariesMaxLevel::refineToMaxLevel(GridStorage& storage,
    RefinementFunctor& functor, unsigned int maxLevel) {
  if (storage.getSize() == 0) {
    throw generation_exception("storage empty");
  }

  // Algorithm should be able to look for several points in grid to refine
  // So we store an array with refinements_num maximal points
  size_t refinements_num = functor.getRefinementsNum();
  RefinementFunctor::value_type* max_values = new
  RefinementFunctor::value_type[refinements_num];
  size_t* max_indexes = new size_t[refinements_num];

  for (size_t i = 0; i < refinements_num; i++) {
    max_values[i] = functor.start();
    max_indexes[i] = 0;
  }

  size_t min_idx = 0;

  RefinementFunctor::value_type max_value = max_values[min_idx];
  size_t max_index = max_indexes[min_idx];

  GridPoint point;
  GridStorage::grid_map_iterator end_iter = storage.end();

  // I think this may be dependent on local support
  for (GridStorage::grid_map_iterator iter = storage.begin(); iter != end_iter;
       iter++) {
    point = *(iter->first);

    GridStorage::grid_map_iterator child_iter;
    bool refineCandidate = false;

    // check if point is on max level in every dimension
    for (size_t dTest = 0; dTest < storage.getDimension(); dTest++) {
      index_t source_index;
      level_t source_level;
      point.get(dTest, source_level, source_index);

      if (source_level < maxLevel) {
        refineCandidate = true;
      }
    }

    // DEBUG
    // std::cout << index.toString() << " " << refineCandidate << std::endl;

    if (refineCandidate == true) {
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
            RefinementFunctor::value_type current_value =
              functor(storage, iter->second);

            // DEBUG
            // std::cout << "iter-second: " << iter->second <<
            // " current_value: " << current_value << std::endl;
            if (current_value > max_value) {
              // Replace the minimal point in result array,
              // find the new minimal point
              max_values[min_idx] = current_value;
              max_indexes[min_idx] = iter->second;
              min_idx = getIndexOfMin(max_values, refinements_num);
              max_value = max_values[min_idx];
              break;
            }
          }
        } else {
          // left child
          point.set(d, source_level + 1, 2 * source_index - 1);
          child_iter = storage.find(&point);

          // if there no more grid points --> test if we should refine the grid
          if (child_iter == end_iter) {
            RefinementFunctor::value_type current_value =
              functor(storage, iter->second);

            // DEBUG
            // std::cout << "iter-second: " << iter->second <<
            // " current_value: " << current_value << std::endl;
            if (current_value > max_value) {
              // Replace the minimal point in result array,
              // find the new  minimal point
              max_values[min_idx] = current_value;
              max_indexes[min_idx] = iter->second;
              min_idx = getIndexOfMin(max_values, refinements_num);
              max_value = max_values[min_idx];
              break;
            }
          }

          // right child
          point.set(d, source_level + 1, 2 * source_index + 1);
          child_iter = storage.find(&point);

          if (child_iter == end_iter) {
            RefinementFunctor::value_type current_value =
              functor(storage, iter->second);

            // DEBUG
            // std::cout << "iter-second: " << iter->second <<
            // " current_value: " << current_value << std::endl;
            if (current_value > max_value) {
              // Replace the minimal point in result array,
              // find the new minimal point
              max_values[min_idx] = current_value;
              max_indexes[min_idx] = iter->second;
              min_idx = getIndexOfMin(max_values, refinements_num);
              max_value = max_values[min_idx];
              break;
            }
          }
        }

        point.set(d, source_level, source_index);
      }
    }
  }

  // DEBUG
  // std::cout << "Num refinements: "  << refinements_num << std::endl;

  // can refine grid on several points
  double threshold = functor.getRefinementThreshold();

  for (size_t i = 0; i < refinements_num; i++) {
    max_value = max_values[i];
    max_index = max_indexes[i];

    // DEBUG
    // std::cout << "Num: " << i << " Max-value: " << max_value << std::endl;
    if (max_value != functor.start() && fabs(max_value) >= threshold) {
      // DEBUG
      // std::cout << "Start refining..." << std::endl;
      refineGridpoint(storage, max_index, maxLevel);
    }
  }

  delete[] max_values;
  delete[] max_indexes;
}


size_t HashRefinementBoundariesMaxLevel::getNumberOfRefinablePointsToMaxLevel(
  GridStorage& storage, unsigned int maxLevel) {
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
    bool refineCandidate = false;

    // check if point is on max level in every dimension
    for (size_t dTest = 0; dTest < storage.getDimension(); dTest++) {
      index_t source_index;
      level_t source_level;
      point.get(dTest, source_level, source_index);

      if (source_level < maxLevel) {
        refineCandidate = true;
      }
    }

    if (refineCandidate == true) {
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
  }

  // DEBUG
  // std::cout << counter << std::endl;

  return counter;
}


void HashRefinementBoundariesMaxLevel::refineGridpoint1D(GridStorage& storage,
                                                         GridPoint& point, size_t d,
                                                         unsigned int maxLevel) {
  index_t source_index;
  level_t source_level;
  point.get(d, source_level, source_index);

  if (source_level < maxLevel) {
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
}


void HashRefinementBoundariesMaxLevel::refineGridpoint(GridStorage& storage,
    size_t refine_index, unsigned int maxLevel) {
  GridPoint point(storage[refine_index]);

  // Sets leaf property of point, which is refined to false
  storage[refine_index].setLeaf(false);

  for (size_t d = 0; d < storage.getDimension(); d++) {
    refineGridpoint1D(storage, point, d, maxLevel);
  }
}


}  // namespace base
}  // namespace sgpp
