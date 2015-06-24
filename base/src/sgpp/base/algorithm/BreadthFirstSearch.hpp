// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef BREADTHFIRSTSEARCH_HPP
#define BREADTHFIRSTSEARCH_HPP

#include <queue>

#include <sgpp/globaldef.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>

namespace SGPP {
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
         * @param storage   storage of the grid
         */
        BreadthFirstSearch(GridStorage* storage) :
          functor(),
          storage(storage) {
        }

        /**
         * Constructor.
         *
         * @param functor   functor to apply in the BFS
         * @param storage   storage of the grid
         */
        BreadthFirstSearch(FUNC& functor, GridStorage* storage) :
          functor(functor),
          storage(storage) {
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
          const size_t n = storage->size();
          const size_t d = storage->dim();

          std::vector<bool> visited(n, false);
          grid_iterator iterator(storage);
          std::queue<size_t> queue;
          size_t index = iterator.seq();
          queue.push(index);
          visited[index] = true;

          while (!queue.empty()) {
            index = queue.front();
            queue.pop();

            const GridIndex& point = *(*storage)[index];
            iterator.set(point);
            functor(source, result, iterator);

            for (size_t t = 0; t < d; t++) {
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

              iterator.up(t);
            }
          }
        }

      protected:
        /// functor to apply in the BFS
        FUNC functor;
        /// storage of the grid
        GridStorage* storage;
    };

  }
}

#endif /* BREADTHFIRSTSEARCH_HPP */
