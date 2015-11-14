// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef HASHGRIDSTORAGE_HPP
#define HASHGRIDSTORAGE_HPP

#include <unordered_map>

#include <memory>
#include <vector>
#include <string>
#include <sstream>
#include <exception>
#include <list>
#include <typeinfo>
#include <stdint.h>

#include <sgpp/globaldef.hpp>


#include <sgpp/base/exception/generation_exception.hpp>

#include <sgpp/base/grid/storage/hashmap/HashGridIndex.hpp>
#include <sgpp/base/grid/storage/hashmap/SerializationVersion.hpp>

#include <sgpp/base/grid/common/BoundingBox.hpp>
#include <sgpp/base/grid/common/Stretching.hpp>

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataMatrixSP.hpp>




namespace SGPP {
  namespace base {

    class HashGridIterator;

    /**
     * Generic hash table based index storage.
     */
    class HashGridStorage {
      public:
        /// type of grid points
        typedef HashGridIndex index_type;
        /// pointer to index_type
        typedef HashGridIndex* index_pointer;
        /// pointer to constant index_type
        typedef const HashGridIndex* index_const_pointer;
        /// unordered_map of index_pointers
        typedef std::unordered_map<index_pointer, size_t, HashGridIndexPointerHashFunctor, HashGridIndexPointerEqualityFunctor > grid_map;
        /// iterator of grid_map
        typedef grid_map::iterator grid_map_iterator;
        /// const_iterator of grid_map
        typedef grid_map::const_iterator grid_map_const_iterator;

        /// vector of index_pointers
        typedef std::vector<index_pointer> grid_list;
        /// iterator of grid_list
        typedef grid_list::iterator grid_list_iterator;

        /// iterator for grid points
        typedef HashGridIterator grid_iterator;

        /**
         * Constructor
         *
         * initializes the boundingBox with a trivial cube
         *
         * @param dim the dimension of the sparse grid
         */
        HashGridStorage(size_t dim);

        /**
         * Constructor
         *
         * initializes the boundingBox with a reference to another boundingbox
         *
         * @param creationBoundingBox reference to bounding box object that describes the grid boundaries
         */
        HashGridStorage(BoundingBox& creationBoundingBox);

        /**
         * Constructor
         *
         * initializes the stretching with a reference to another stretching
         *
         * @param creationStretching reference to stretching object that describes the grid boundaries and the stretching
         */
        HashGridStorage(Stretching& creationStretching);

        /**
         * Constructor that reads the data from a string
         *
         * @param istr the string that contains the data
         */
        HashGridStorage(std::string& istr);

        /**
         * Constructor that reads the data from an input stream
         *
         * @param istream the inputstream that contains the data
         */
        HashGridStorage(std::istream& istream);

        /**
         * Copy Constructor
         */
        HashGridStorage(HashGridStorage& copyFrom);

        /**
         * Destructor
         */
        ~HashGridStorage();

        /**
         * deletes all grid points in the storage
         */
        void emptyStorage();

        /**
         * Remove several point from HashGridStorage. The points to removed
         * are stored in a list. This function returns a vector of remaining points
         * given by their
         *
         * @param removePoints vector containing the indices of the points that should be removed
         *
         * @return a vector containing the indices of remaining points given by their "old" index
         */
        std::vector<size_t> deletePoints(std::list<size_t>& removePoints);

        /**
         * unserializes the grid from a string, algorithmic dimensions are not reseted
         *
         * @param istr the string that contains the data
         */
        void unserialize_noAlgoDims(std::string& istr);

        /**
         * serialize the gridstorage into a string
         *
         * @return a string that contains all gridstorage information
         */
        std::string serialize();

        /**
         * serialize the gridstorage into a stream
         *
         * @param ostream reference to a stream into that all gridstorage information is written
         */
        void serialize(std::ostream& ostream);

        /**
         * serialize the gridstorage's gridpoints into a stream
         *
         * @return returns the string that contains all gridpoint information
         */
        std::string toString();

        /**
         * serialize the gridstorage's gridpoints into a stream
         *
         * @param stream reference to a stream into that all gridpoint information is written
         */
        void toString(std::ostream& stream);

        /**
         * gets the size of the hashmap
         *
         * @return returns the size of the hashmap
         */
        size_t size() const;

        /**
         * gets the number of inner grid points
         *
         * @return the number of inner grid points
         */
        size_t getNumInnerPoints() const;

        /**
         * gets the dimension of the grid
         *
         * @return the dimension of the grid stored in this HashGridStorage object
         */
        size_t dim() const;

        /**
         * gets the index number for given gridpoint by its sequence number
         *
         * @param seq the sequence number of the index
         *
         * @return gridindex object (pointer)
         */
        inline index_pointer operator[](size_t seq) {
          return list[seq];
        }

        /**
         * gets the index number for given gridpoint by its sequence number
         *
         * @param seq the sequence number of the index
         * @return gridindex object (constant pointer)
         */
        inline index_const_pointer operator[](size_t seq) const {
          return list[seq];
        }

        /**
         * gets the index number for given gridpoint by its sequence number
         *
         * @param seq the sequence number of the index
         *
         * @return gridindex object (pointer)
         */
        inline index_pointer get(size_t seq) const {
          return list[seq];
        }

        /**
         * insert a new index into map
         *
         * @param index reference to the index that should be inserted
         *
         * @return
         */
        size_t insert(index_type& index);

        /**
         * updates an already stored index
         *
         * @param index reference to the index that should be updated
         * @param pos position where the index should be stored
         */
        void update(index_type& index, size_t pos);

        /**
         * This methods removes the gridpoint added last. Use with coution, only needed for
         * expanding the grid because of the shadow-storage of prewavelets. Please refer to the
         * Prewavelet grid for further description of the shadow storage.
         *
         */
        void deleteLast();

        /**
         * creates a pointer to index from a reference to index by creating
         * a new instance of a index object
         *
         * @param index address of index object
         *
         * @return pointer to new index object
         */
        index_pointer create(index_type& index);

        /**
         * removes an index from gridstorage
         *
         * @param index pointer to index that should be removed
         */
        void destroy(index_pointer index);

        /**
         * stores a given index in the hashmap
         *
         * @param index pointer to index that should be stored
         *
         * @return sequence number
         */
        unsigned int store(index_pointer index);

        /**
         * sets the iterator to a given index
         *
         * @param index the index to which the cursor should be moved
         * @return iterator pointing to the index
         */
        grid_map_iterator find(index_pointer index);

        /**
         * set iterator to the first position in the map
         * @return iterator pointing to the beginning of the map
         */
        grid_map_iterator begin();

        /**
         * sets the iterator to last position in the map
         * @return iterator pointing to the end of the map
         */
        grid_map_iterator end();

        /**
         * Tests if index is in the storage
         *
         * @param index pointer to index that should be tested
         *
         * @return true if the index is in the storage
         */
        bool has_key(HashGridIndex* index);

        /**
         * Sets the index to the left level zero parent
         *
         * @param index pointer to index the should be modified
         * @param dim the dimension in which the modification is taken place
         */
        void left_levelzero(HashGridIndex* index, size_t dim);

        /**
         * Sets the index to the right level zero parent
         *
         * @param index pointer to index the should be modified
         * @param dim the dimension in which the modification is taken place
         */
        void right_levelzero(HashGridIndex* index, size_t dim);

        /**
         * Sets the index to the left child
         *
         * @param index pointer to index the should be modified
         * @param dim the dimension in which the modification is taken place
         */
        void left_child(HashGridIndex* index, size_t dim);

        /**
         * Sets the index to the right child
         *
         * @param index pointer to index the should be modified
         * @param dim the dimension in which the modification is taken place
         */
        void right_child(HashGridIndex* index, size_t dim);

        /**
         * Resets the index to the top level in direction d
         *
         * @param index pointer to index the should be modified
         * @param d the dimension in which the modification is taken place
         */
        void top(HashGridIndex* index, size_t d);

        /**
         * Gets the seq number for index
         *
         * @param index pointer to index which sequence number should be determined
         *
         * @return the seq number for index
         */
        size_t seq(HashGridIndex* index);

        /**
         * Tests if seq number does not point to a valid grid index
         *
         * @param s sequence number that should be tested
         *
         * @return true if we are not EOF
         */
        bool end(size_t s);

        /**
         * returns the algorithmic dimensions (the dimensions in which the Up Down
         * operations should be applied)
         *
         * @return the algorithmic dimensions
         */
        std::vector<size_t> getAlgorithmicDimensions();

        /**
         * sets the algorithmic dimensions (the dimensions in which the Up Down
         * operations should be applied)
         *
         * @param newAlgoDims std::vector containing the algorithmic dimensions
         */
        void setAlgorithmicDimensions(std::vector<size_t> newAlgoDims);

        /**
         * Recalculates the leaf-property of every grid point.
         * This might be useful in case of a grid unserialization
         */
        void recalcLeafProperty();

        /**
         * get the bounding box of the current grid
         *
         * @return returns a pointer to HashGridStorage's bounding box
         */
        BoundingBox* getBoundingBox();

        /**
         * get the stretching bounding box of the current grid
         *
         * @return returns a pointer to HashGridStorage's bounding box
         */
        Stretching* getStretching();

        /**
         * sets the bounding box of the current grid
         *
         * @param bb bounding box to which the HashGridStorage's pointer is set
         */
        void setBoundingBox(BoundingBox& bb);

        /**
         * sets the stretching bounding box of the current grid
         *
         * @param bb stretching to which the HashGridStorage's pointer is set
         */
        void setStretching(Stretching& bb);

        /**
         * Converts this storage from AOS (array of structures) to SOA (structure of array)
         * with modification to speed up iterative function evaluation. The Level
         * array won't contain the levels, it contains the level to the power of two
         *
         * @param level DataMatrix to store the grid's level to the power of two
         * @param index DataMatrix to store the grid's indices
         */
        void getLevelIndexArraysForEval(DataMatrix& level, DataMatrix& index);

        /**
         * Converts this storage from AOS (array of structures) to SOA (structure of array)
         * with modification to speed up iterative function evaluation. The Level
         * array won't contain the levels, it contains the level to the power of two
         *
         * @param level DataMatrix to store the grid's level to the power of two
         * @param index DataMatrix to store the grid's indices
         */
        void getLevelIndexArraysForEval(DataMatrixSP& level, DataMatrixSP& index);

        /**
         * Converts this storage from AOS (array of structures) to SOA (structure of array)
         * with modification to speed up iterative Laplace Calculations: the level
         * won't contain the levels, it contains 2 to the neagative power of the level.
         *
         * @param level DataMatrix to store the grid's modified level
         */
        void getLevelForIntegral(DataMatrix& level);


        /**
         * returns the max. depth in all dimension of the grid
         */
        size_t getMaxLevel() const;

        /**
         * Converts this storage from AOS (array of structures) to SOA (structure of array)
         * with modification to speed up iterative function evaluation. The Level
         * array won't contain the levels, it contains the level to the power of two.
         *
         * The returned format is only useful for a multi-evaluation of modlinear grids
         *
         * @param level DataMatrix to store the grid's level to the power of two
         * @param index DataMatrix to store the grid's indices
         * @param mask DataMatrix to store masks of operations
         * @param offset DataMatrix to store offset for operations
         */
        void getLevelIndexMaskArraysForModEval(DataMatrix& level, DataMatrix& index,
                                               DataMatrix& mask, DataMatrix& offset);

        /**
         * Converts this storage from AOS (array of structures) to SOA (structure of array)
         * with modification to speed up iterative function evaluation. The Level
         * array won't contain the levels, it contains the level to the power of two.
         *
         * The returned format is only useful for a multi-evaluation of modlinear grids
         *
         * @param level DataMatrixSP to store the grid's level to the power of two
         * @param index DataMatrixSP to store the grid's indices
         * @param mask DataMatrixSP to store masks of operations
         * @param offset DataMatrixSP to store offset for operations
         */
        void getLevelIndexMaskArraysForModEval(DataMatrixSP& level, DataMatrixSP& index,
                                               DataMatrixSP& mask, DataMatrixSP& offset);

      protected:
        /**
         * returns the next sequence numbers
         *
         * @return returns the next sequence numbers
         */
        size_t seq() const;

      private:

        /// the dimension of the grid
        size_t DIM;

        /// the grid points
        grid_list list;
        /// the indices of the grid points
        grid_map map;
        /// algorithmic dimension, these are used in Up/Downs
        std::vector<size_t> algoDims;

        /// the grid's bounding box
        BoundingBox* boundingBox;
        /// the grid's stretching
        Stretching* stretching;

        /// Flag to check if stretching or bb used
        bool bUseStretching;



        /**
         * Parses the gird's information (grid points, dimensions, bounding box) from a string stream
         *
         * @param istream the string stream that contains the information
         */
        void parseGridDescription(std::istream& istream);

    };




    HashGridStorage::index_pointer
    inline HashGridStorage::create(index_type& index) {
      index_pointer insert = new HashGridIndex(index);
      return insert;
    }

    void
    inline HashGridStorage::destroy(index_pointer index) {
      delete index;
    }

    unsigned int
    inline HashGridStorage::store(index_pointer index) {
      list.push_back(index);
      return static_cast<unsigned int>(map[index] =
                                         static_cast<unsigned int>(this->seq() - 1));
    }

    HashGridStorage::grid_map_iterator
    inline HashGridStorage::find(index_pointer index) {
      return map.find(index);
    }

    HashGridStorage::grid_map_iterator
    inline HashGridStorage::begin() {
      return map.begin();
    }

    HashGridStorage::grid_map_iterator
    inline HashGridStorage::end() {
      return map.end();
    }

    bool
    inline HashGridStorage::has_key(HashGridIndex* index) {
      return map.find(index) != map.end();
    }

    void
    inline HashGridStorage::left_levelzero(HashGridIndex* index, size_t dim) {
      index_type::level_type l;
      index_type::index_type i;
      index->get(dim, l, i);
      index->set(dim, 0, 0);
    }

    void
    inline HashGridStorage::right_levelzero(HashGridIndex* index, size_t dim) {
      index_type::level_type l;
      index_type::index_type i;
      index->get(dim, l, i);
      index->set(dim, 0, 1);
    }

    void
    inline HashGridStorage::left_child(HashGridIndex* index, size_t dim) {
      index_type::level_type l;
      index_type::index_type i;
      index->get(dim, l, i);
      index->set(dim, l + 1, 2 * i - 1);
    }

    void
    inline HashGridStorage::right_child(HashGridIndex* index, size_t dim) {
      index_type::level_type l;
      index_type::index_type i;
      index->get(dim, l, i);
      index->set(dim, l + 1, 2 * i + 1);
    }

    void
    inline HashGridStorage::top(HashGridIndex* index, size_t d) {
      index->set(d, 1, 1);
    }

    size_t
    inline HashGridStorage::seq(HashGridIndex* index) {
      grid_map_iterator iter = map.find(index);

      if (iter != map.end()) {
        return iter->second;
      } else {
        return map.size() + 1;
      }
    }

    bool
    inline HashGridStorage::end(size_t s) {
      return s > map.size();
    }

    std::vector<size_t>
    inline HashGridStorage::getAlgorithmicDimensions() {
      return algoDims;
    }



    size_t
    inline HashGridStorage::seq() const {
      return list.size();
    }




  } // namespace base
} // namespace SGPP

#include <sgpp/base/grid/storage/hashmap/HashGridIterator.hpp>


#endif /* HASHGRIDSTORAGE_HPP */
