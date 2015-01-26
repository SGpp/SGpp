/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de), Dirk Pflueger (pflueged@in.tum.de)

#ifndef HASHGRIDITERATOR_HPP
#define HASHGRIDITERATOR_HPP

#include "base/tools/hash_map_config.hpp"

#include "base/exception/generation_exception.hpp"

#include "base/grid/storage/hashmap/HashGridIndex.hpp"
#include "base/grid/GridStorage.hpp"
#include "base/grid/storage/hashmap/SerializationVersion.hpp"

#include <memory>
#include <string>
#include <sstream>
#include <exception>

namespace sg {
  namespace base {

    template<typename GIT>
    class HashGridStorage;

    /**
     * This class can be used for storage agnostic algorithms.
     * GridIndex has to support: constructor, get, set, push, rehash
     */
    template<typename GIT>
    class HashGridIterator {
      public:
        typedef GIT index_type;
        typedef typename GIT::index_type index_t;
        typedef typename GIT::level_type level_t;

        /**
         * Constructor of the griditerator object
         *
         * @param storage pointer the hashmap that stores the grid points
         */
        HashGridIterator(HashGridStorage<GIT>* storage) :
          storage(storage), index(storage->dim()) {
          for (size_t i = 0; i < storage->dim(); i++) {
            index.push(i, 1, 1);
          }

          index.rehash();
          this->seq_ = storage->seq(&index);
        }

        /**
         * Copy Constructor of the griditerator object
         *
         * @param copy a HashGridIterator object that is used to build this instance
         */
        HashGridIterator(HashGridIterator<GIT>& copy) :
          storage(copy.storage), index(copy.storage->dim()) {
          typename index_type::level_type l;
          typename index_type::index_type i;

          for (size_t dim = 0; dim < storage->dim(); dim++) {
            copy.get(dim, l, i);
            index.push(dim, l, i);
          }

          index.rehash();
          this->seq_ = storage->seq(&index);
        }

        /**
         *  Sets 0,0 in every dimension (Left Level zero ansatzfunction)
         */
        void resetToLevelZero() {
          for (size_t i = 0; i < storage->dim(); i++) {
            index.push(i, 0, 0);
          }

          index.rehash();
          this->seq_ = storage->seq(&index);
        }

        /**
         * left level zero ansatz function for a given dimension
         *
         * @param dim dimension in which we should step to level zero
         */
        void left_levelzero(size_t dim) {
          index.set(dim, 0, 0);
          this->seq_ = storage->seq(&index);
        }

        /**
         * right level zero ansatz function for a given dimension
         *
         * @param dim dimension in which we should step to level zero
         */
        void right_levelzero(size_t dim) {
          index.set(dim, 0, 1);
          this->seq_ = storage->seq(&index);
        }

        /**
         * left child in direction dim
         *
         * @param dim dimension in which we should step to the left child
         */
        void left_child(size_t dim) {
          typename index_type::level_type l;
          typename index_type::index_type i;
          index.get(dim, l, i);
          index.set(dim, l + 1, 2 * i - 1);
          this->seq_ = storage->seq(&index);
        }

        /**
         * right child in direction dim
         *
         * @param dim dimension in which we should step to the right child
         */
        void right_child(size_t dim) {
          typename index_type::level_type l;
          typename index_type::index_type i;
          index.get(dim, l, i);
          index.set(dim, l + 1, 2 * i + 1);
          this->seq_ = storage->seq(&index);
        }

        /**
         * resets the iterator to the top if dimension d
         *
         * @todo (heinecke, must) maybe rename to steptoLevelOne
         *
         * @param d the moving direction
         */
        void top(size_t d) {
          index.set(d, 1, 1);
          this->seq_ = storage->seq(&index);
        }

        /**
         * hierarchical parent in direction dim
         *
         * @param d the moving direction
         */
        void up(size_t d) {
          typename index_type::level_type l;
          typename index_type::index_type i;
          index.get(d, l, i);

          i /= 2;
          i += i % 2 == 0 ? 1 : 0;

          index.set(d, l - 1, i);
          this->seq_ = storage->seq(&index);
        }

        /**
         * step left in direction dim
         *
         * @param d the moving direction
         */
        void step_left(size_t d) {
          typename index_type::level_type l;
          typename index_type::index_type i;
          index.get(d, l, i);
          index.set(d, l, i - 2);
          this->seq_ = storage->seq(&index);

        }

        /**
         * step right in direction dim
         *
         * @param d the moving direction
         */
        void step_right(size_t d) {
          typename index_type::level_type l;
          typename index_type::index_type i;
          index.get(d, l, i);
          index.set(d, l, i + 2);
          this->seq_ = storage->seq(&index);

        }

        bool isInnerPoint() {
          return index.isInnerPoint();
        }

        /**
         * returns true if there are no more childs in any dimension
         *
         * @return returns true if there are no more childs in any dimension
         */
        bool hint() const {
          return storage->get(this->seq_)->isLeaf();
        }

        /**
         * returns true if there are more left childs in dimension d
         *
         * @param d the moving direction
         */
        bool hint_left(size_t d) {
          typename index_type::level_type l;
          typename index_type::index_type i;
          bool hasIndex = true;

          index.get(d, l, i);
          index.set(d, l + 1, 2 * i - 1);

          GIT* my_Index = index.getPointer();
          hasIndex = storage->has_key(my_Index);

          index.set(d, l, i);

          return hasIndex;
        }

        /**
         * returns true if there are more right childs in dimension d
         *
         * @param d the moving direction
         */
        bool hint_right(size_t d) {
          typename index_type::level_type l;
          typename index_type::index_type i;
          bool hasIndex = true;

          index.get(d, l, i);
          index.set(d, l + 1, 2 * i + 1);

          GIT* my_Index = index.getPointer();
          hasIndex = storage->has_key(my_Index);

          index.set(d, l, i);

          return hasIndex;
        }

        /**
         * Gets level @c l and index @c i in dimension @c d of the current grid point
         *
         * @param d the dimension of interest
         * @param l the ansatz function's level
         * @param i the ansatz function's index
         */
        void get(size_t d, typename index_type::level_type& l,
                 typename index_type::index_type& i) const {
          index.get(d, l, i);
        }

        /**
         * Sets level @c l and index @c i in dimension @c d of the current grid point.
               * Recomputes the hash value of the current grid point.
         *
         * @param d the dimension of interest
         * @param l the ansatz function's level
         * @param i the ansatz function's index
         */
        void set(size_t d, typename index_type::level_type l,
                 typename index_type::index_type i) {
          index.set(d, l, i);
          this->seq_ = storage->seq(&index);
        }

        /**
         * Sets level @c l and index @c i in dimension @c d of the current grid point.
               * Does not recompute hash value of the current grid point.
         *
         * @param d the dimension of the gridpoint
         * @param l the ansatz function's level
         * @param i the ansatz function's index
         */
        void push(size_t d, typename index_type::level_type l,
                  typename index_type::index_type i) {
          index.push(d, l, i);
        }

        /**
         * returns the current sequence number
         *
         * @return the current sequence number
         */
        size_t seq() const {
          return seq_;
        }

        /**
         * Returns the the maximal level of the grid in the given dimension.
         *
         * @param dim the dimension
         */
        level_t getGridDepth(size_t dim) {

          typename index_type::level_type depth = 1;
          typename index_type::level_type orig_level, cur_level;
          typename index_type::index_type orig_index, cur_index;

          index.get(dim, orig_level, orig_index);

          while (true) {
            if (this->hint_left(dim)) {
              depth++;
              this->left_child(dim);
            } else if (this->hint_right(dim)) {
              depth++;
              this->right_child(dim);
            } else {

              index.get(dim, cur_level, cur_index);

              bool hasFound = false; //Was a next index found?

              //Ok, we have no more childs left. Now we slide from left to right in the dim on
              //the same level, to see, if there are adaptive refinements
              for (size_t i = cur_index + 2; i
                   < (unsigned int) (1 << (depth)); i = i + 2) {
                this->set(dim, cur_level, *reinterpret_cast< unsigned int* >(&i));

                //does this index exist?
                if (!storage->end(this->seq())) {
                  if (this->hint_left(dim)) {
                    depth++;
                    this->left_child(dim);
                    hasFound = true;
                    break;
                  } else if (this->hint_right(dim)) {
                    depth++;
                    this->right_child(dim);
                    hasFound = true;
                    break;
                  }
                }
              }

              if (!hasFound) {
                break;
              }

            }
          }

          this->set(dim, orig_level, orig_index);
          return depth;
        }

        std::string toString() {
          return index.toString();
        }

      private:
        /// Pointer the the hashmap that stores the gridpoints
        HashGridStorage<GIT>* storage;
        /// GridIndex object used to operate on the current position in the hashmap
        GIT index;
        /// true if the current point is a leaf, otherwise false
        //bool Leaf;
        /// the current gridpoint's index
        size_t seq_;
    };

  }
}

#endif /* HASHGRIDITERATOR_HPP */
