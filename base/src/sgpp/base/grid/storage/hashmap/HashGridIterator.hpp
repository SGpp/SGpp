// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef HASHGRIDITERATOR_HPP
#define HASHGRIDITERATOR_HPP

#include <sgpp/base/exception/generation_exception.hpp>

#include <sgpp/base/grid/storage/hashmap/HashGridIndex.hpp>
#include <sgpp/base/grid/storage/hashmap/HashGridStorage.hpp>

#include <memory>
#include <string>
#include <sstream>
#include <exception>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * This class can be used for storage agnostic algorithms.
     * GridIndex has to support: constructor, get, set, push, rehash
     */
    class HashGridIterator {
      public:
        /// type of grid points
        typedef HashGridIndex index_type;
        /// index type
        typedef HashGridIndex::index_type index_t;
        /// level type
        typedef HashGridIndex::level_type level_t;

        /**
         * Constructor of the griditerator object
         *
         * @param storage pointer the hashmap that stores the grid points
         */
        HashGridIterator(HashGridStorage* storage);

        /**
         * Copy Constructor of the griditerator object
         *
         * @param copy a HashGridIterator object that is used to build this instance
         */
        HashGridIterator(HashGridIterator& copy);

        /**
         * Destructor
         */
        ~HashGridIterator();

        /**
         *  Sets 0,0 in every dimension (Left Level zero ansatzfunction)
         */
        void resetToLevelZero();

        /**
         * left level zero ansatz function for a given dimension
         *
         * @param dim dimension in which we should step to level zero
         */
        void resetToLeftLevelZero(size_t dim);

        /**
         * right level zero ansatz function for a given dimension
         *
         * @param dim dimension in which we should step to level zero
         */
        void resetToRightLevelZero(size_t dim);

        /**
         * resets the iterator to the top if dimension d
         *
         * @param d the moving direction
         */
        void resetToLevelOne(size_t d);

        /**
         * left child in direction dim
         *
         * @param dim dimension in which we should step to the left child
         */
        void leftChild(size_t dim);

        /**
         * right child in direction dim
         *
         * @param dim dimension in which we should step to the right child
         */
        void rightChild(size_t dim);

        /**
         * hierarchical parent in direction dim
         *
         * @param d the moving direction
         */
        void up(size_t d);

        /**
         * step left in direction dim
         *
         * @param d the moving direction
         */
        void stepLeft(size_t d);

        /**
         * step right in direction dim
         *
         * @param d the moving direction
         */
        void stepRight(size_t d);

        /**
         * determines if the grid point is an inner grid point
         *
         * @return true if the grid point is an inner grid point
         */
        bool isInnerPoint() const;

        /**
         * returns true if there are no more children in any dimension
         *
         * @return returns true if there are no more children in any dimension
         */
        bool hint() const;

        /**
         * returns true if there are more left children in dimension d
         *
         * @param d the moving direction
         * @return true if there are more left children in dimension d
         */
        bool hintLeft(size_t d);

        /**
         * returns true if there are more right children in dimension d
         *
         * @param d the moving direction
         * @return true if there are more right children in dimension d
         */
        bool hintRight(size_t d);

        /**
         * Gets level @c l and index @c i in dimension @c d of the current grid point
         *
         * @param d the dimension of interest
         * @param l the ansatz function's level
         * @param i the ansatz function's index
         */
        void get(size_t d, index_type::level_type& l, index_type::index_type& i) const;

        /**
         * Sets level @c l and index @c i in dimension @c d of the current grid point.
         * Recomputes the hash value of the current grid point.
         *
         * @param d the dimension of interest
         * @param l the ansatz function's level
         * @param i the ansatz function's index
         */
        void set(size_t d, index_type::level_type l,
                 index_type::index_type i);

        /**
         * Sets the iterator to a grid point.
         * Recomputes the hash value of the current grid point.
         *
         * @param point new grid point
         */
        void set(const index_type& point);

        /**
         * Sets level @c l and index @c i in dimension @c d of the current grid point.
         * Does not recompute hash value of the current grid point.
         *
         * @param d the dimension of the gridpoint
         * @param l the ansatz function's level
         * @param i the ansatz function's index
         */
        void push(size_t d, index_type::level_type l, index_type::index_type i);

        /**
         * returns the current sequence number
         *
         * @return the current sequence number
         */
        size_t seq() const;

        /**
         * Returns the the maximal level of the grid in the given dimension.
         *
         * @param dim the dimension
         */
        level_t getGridDepth(size_t dim);

        /**
         * Generates a string with level and index of the gridpoint.
         *
         * @returns string into which the gridpoint is written
         */
        std::string toString();

      private:
        /// Pointer the the hashmap that stores the gridpoints
        HashGridStorage* storage;
        /// GridIndex object used to operate on the current position in the hashmap
        HashGridIndex index;
        /// true if the current point is a leaf, otherwise false
        bool Leaf;
        /// the current gridpoint's index
        size_t seq_;
    };

  } // namespace base
} // namespace SGPP

#endif /* HASHGRIDITERATOR_HPP */
