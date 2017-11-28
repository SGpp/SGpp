// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef BREADTHFIRSTSEARCH_HPP
#define BREADTHFIRSTSEARCH_HPP

#include <sgpp/globaldef.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>

#include <queue>
#include <vector>

namespace sgpp {
namespace base {

/**
 * Class for applying a functor of type FUNC on all grid points of a
 * sparse grid in a breadth-first search (BFS) manner.
 * Especially useful for hierarchization with basis functions which
 * fulfill a Lagrange property, e.g., fundamental splines.
 * The functor should implement
 * void operator()(const DataVector& source,
 *                 DataVector& result,
 *                 const grid_iterator& iterator),
 * where source is the original data (e.g., node values),
 * result is the result of the functor (e.g., surplusses), and
 * iterator is at the current grid point.
 */
template<class FUNC>
class BreadthFirstSearch {
 protected:
  /// grid iterator
  typedef GridStorage::grid_iterator grid_iterator;

 public:
  /**
   * Constructor.
   * Uses the default constructor for the functor.
   *
   * @param storage           storage of the grid
   * @param startAtCorners    whether to start the BFS at the corners instead at the midpoint
   */
  explicit BreadthFirstSearch(GridStorage& storage, bool startAtCorners) :
    functor(),
    storage(storage),
    startAtCorners(startAtCorners) {
  }

  /**
   * Constructor.
   *
   * @param functor           functor to apply in the BFS
   * @param storage           storage of the grid
   * @param startAtCorners    whether to start the BFS at the corners instead at the midpoint
   */
  BreadthFirstSearch(FUNC& functor, GridStorage& storage, bool startAtCorners) :
    functor(functor),
    storage(storage),
    startAtCorners(startAtCorners) {
  }

  /**
   * Destructor.
   */
  ~BreadthFirstSearch() {
  }

  /**
   * Execute the BFS.
   *
   * @param[in]  source   original data (e.g., node values)
   * @param[out] result   result of the BFS (e.g., surplusses)
   */
  template<class DataType>
  void execute(const DataType& source, DataType& result) {
    const size_t n = storage.getSize();
    const size_t d = storage.getDimension();

    std::vector<bool> visited(n, false);
    grid_iterator iterator(storage);
    std::queue<size_t> queue;
    size_t index;
    level_t l;
    index_t i;

    for (size_t q = 0; q < 2; q++) {
      if (startAtCorners && (q == 0)) {
        // insert all corners of [0, 1]^d into queue
        iterator.resetToLevelZero();

        std::vector<size_t> powersOfTwo(d + 1);
        powersOfTwo[0] = 1;

        for (size_t t = 0; t < d; t++) {
          powersOfTwo[t + 1] = 2 * powersOfTwo[t];
        }

        // iterate through the 2^d corners
        for (size_t k = 0; k < powersOfTwo[d]; k++) {
          for (size_t t = 0; t < d; t++) {
            iterator.set(t, 0, static_cast<index_t>((k & powersOfTwo[t]) != 0));
          }

          // add the corner to the queue if present in the grid
          index = iterator.seq();

          if (index < n) {
            queue.push(index);
            visited[index] = true;
          }
        }
      } else {
        // insert midpoint (0.5, ..., 0.5) into queue
        index = iterator.seq();
        queue.push(index);
        visited[index] = true;
      }

      // work through queue until empty
      while (!queue.empty()) {
        index = queue.front();
        queue.pop();

        if (index >= n) {
          continue;
        }

        const GridPoint& point = storage[index];
        iterator.set(point);
        functor(source, result, iterator);

        for (size_t t = 0; t < d; t++) {
          iterator.get(t, l, i);

          if (l == 0) {
            // boundary point
            iterator.set(t, 1, 1);
            index = iterator.seq();

            if ((index < n) && (!visited[index])) {
              queue.push(index);
              visited[index] = true;
            }
          } else {
            // interior point
            {
              iterator.leftChild(t);
              index = iterator.seq();

              if ((index < n) && (!visited[index])) {
                queue.push(index);
                visited[index] = true;
              }
            }

            {
              iterator.stepRight(t);
              index = iterator.seq();

              if ((index < n) && (!visited[index])) {
                queue.push(index);
                visited[index] = true;
              }
            }
          }

          iterator.set(t, l, i);
        }
      }

      // repeat once iff we started at corners, but did not reach the midpoint
      if (q == 0) {
        if (startAtCorners) {
          for (size_t t = 0; t < d; t++) {
            iterator.resetToLevelOne(t);
          }

          index = iterator.seq();

          if ((index < n) || visited[index]) {
            break;
          }
        } else {
          break;
        }
      }
    }
  }

 protected:
  /// functor to apply in the BFS
  FUNC functor;
  /// storage of the grid
  GridStorage& storage;
  /// whether to start the BFS at the corners instead at the midpoint
  bool startAtCorners;
};

}  // namespace base
}  // namespace sgpp

#endif /* BREADTHFIRSTSEARCH_HPP */
