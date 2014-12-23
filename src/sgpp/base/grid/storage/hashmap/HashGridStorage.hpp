/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Dirk Pflueger (dirk.pflueger@in.tum.de), Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef HASHGRIDSTORAGE_HPP
#define HASHGRIDSTORAGE_HPP

//#include "base/tools/hash_map_config.hpp"

#include <unordered_map>

#include "base/exception/generation_exception.hpp"

#include "base/grid/storage/hashmap/HashGridIndex.hpp"
#include "base/grid/storage/hashmap/HashGridIterator.hpp"
#include "base/grid/storage/hashmap/SerializationVersion.hpp"

#include "base/grid/common/BoundingBox.hpp"
#include "base/grid/common/Stretching.hpp"

#include "base/datatypes/DataMatrix.hpp"
#include "base/datatypes/DataMatrixSP.hpp"

#ifdef SG_PARALLEL
#include "parallel/tools/TypesParallel.hpp"
#include "parallel/datadriven/tools/DMVectorizationPaddingAssistant.hpp"
#endif

#include <memory>
#include <vector>
#include <string>
#include <sstream>
#include <exception>
#include <list>
#include <typeinfo>
#include <stdint.h>

namespace sg {
  namespace base {

    /**
     * Generic hash table based index storage.
     */
    template<typename GIT>
    class HashGridStorage {
      public:
        typedef GIT index_type;
        typedef GIT* index_pointer;
//#ifndef USETRONE
//        typedef std::hash_map<index_pointer, size_t, sg::base::hash<index_pointer>, sg::base::eqIndex<index_pointer> > grid_map;
//#else
        //typedef std::tr1::unordered_map<index_pointer, size_t, hash<index_pointer>, eqIndex<index_pointer> > grid_map;
//#endif
        typedef std::unordered_map<index_pointer, size_t, sg::base::hash<index_pointer>, sg::base::eqIndex<index_pointer> > grid_map;

        typedef typename grid_map::iterator grid_map_iterator;
        typedef typename grid_map::const_iterator grid_map_const_iterator;

        typedef std::vector<index_pointer> grid_list;
        typedef typename grid_list::iterator grid_list_iterator;

        typedef HashGridIterator<GIT> grid_iterator;

        /**
         * Constructor
         *
         * initializes the boundingBox with a trivial cube
         *
         * @param dim the dimension of the sparse grid
         */
        HashGridStorage(size_t dim) : DIM(dim), list(), map() {
          bUseStretching = false;
          boundingBox = new BoundingBox(DIM);
          stretching = NULL;

          for (size_t i = 0; i < DIM; i++) {
            algoDims.push_back(i);
          }
        }

        /**
         * Constructor
         *
         * initializes the boundingBox with a reference to another boundingbox
         *
         * @param creationBoundingBox reference to bounding box object that describes the grid boundaries
         */
        HashGridStorage(BoundingBox& creationBoundingBox) : DIM(), list(), map() {
          bUseStretching = false;
          boundingBox = new BoundingBox(creationBoundingBox);
          stretching = NULL;
          DIM = boundingBox->getDimensions();

          for (size_t i = 0; i < DIM; i++) {
            algoDims.push_back(i);
          }
        }

        /**
         * Constructor
         *
         * initializes the stretching with a reference to another stretching
         *
         * @param creationStretching reference to stretching object that describes the grid boundaries and the stretching
         */
        HashGridStorage(Stretching& creationStretching) : DIM(), list(), map() {
          bUseStretching = true;
          boundingBox = NULL;
          stretching = new Stretching(creationStretching);
          DIM = stretching->getDimensions();

          for (size_t i = 0; i < DIM; i++) {
            algoDims.push_back(i);
          }
        }

        /**
         * Constructor that reads the data from a string
         *
         * @param istr the string that contains the data
         */
        HashGridStorage(std::string& istr) : DIM(0), list(), map() {
          std::istringstream istream;
          istream.str(istr);

          parseGridDescription(istream);

          for (size_t i = 0; i < DIM; i++) {
            algoDims.push_back(i);
          }
        }

        /**
         * Constructor that reads the data from an input stream
         *
         * @param istream the inputstream that contains the data
         */
        HashGridStorage(std::istream& istream) : DIM(0), list(), map() {
          parseGridDescription(istream);

          for (size_t i = 0; i < DIM; i++) {
            algoDims.push_back(i);
          }
        }

        /**
         * Copy Constructor
         */
        HashGridStorage(HashGridStorage& copyFrom) : DIM(0), list(), map() {
          // Copy dimensions and bounding box
          DIM = copyFrom.DIM;

          if (!copyFrom.bUseStretching) {
            bUseStretching = false;
            boundingBox = new BoundingBox(*(copyFrom.boundingBox));
            stretching = NULL;
          }

          if (copyFrom.bUseStretching) {
            bUseStretching = true;
            stretching = new Stretching(*(copyFrom.stretching));
            boundingBox = NULL;
          }

          // copy algorithmic dimensions
          for (size_t i = 0; i > copyFrom.algoDims.size(); i++) {
            algoDims.push_back(copyFrom.algoDims[i]);
          }

          // copy gridpoints
          for (size_t i = 0; i < copyFrom.size(); i++) {
            this->insert(*(copyFrom[i]));
          }
        }


        /**
         * Destructor
         */
        ~HashGridStorage() {
          // delete all grid points
          for (grid_list_iterator iter = list.begin(); iter != list.end(); iter++) {
            delete *iter;
          }

          // delete the grid's bounding box
          if (!bUseStretching) {
            delete boundingBox;
          }

          if (bUseStretching) {
            delete stretching;
          }
        }

        /**
         * deletes all grid points in the storage
         */
        void emptyStorage() {
          // delete all grid points
          for (grid_list_iterator iter = list.begin(); iter != list.end(); iter++) {
            delete *iter;
          }

          // remove all elements from hashmap
          map.clear();
          // remove all list entries
          list.clear();
        }

        /**
         * Remove several point from GridStorage. The points to removed
         * are stored in a list. This function returns a vector of remaining points
         * given by their
         *
         * @param removePoints vector containing the indices of the points that should be removed
         *
         * @return a vector containing the indices of remaining points given by their "old" index
         */
        std::vector<size_t> deletePoints(std::list<size_t>& removePoints) {
          index_pointer curPoint;
          std::vector<size_t> remainingPoints;
          size_t delCounter = 0;

          // sort list
          removePoints.sort();

          //DEBUG : print list points to delete, sorted
          //std::cout << std::endl << "List of points to delete, sorted" << std::endl;
          //for(std::list<size_t>::iterator iter = removePoints.begin(); iter != removePoints.end(); iter++)
          //{
          //  std::cout << " " << *iter << " ";
          //}
          //std::cout << std::endl;

          // Remove points with given indices for index vector and hashmap
          for (std::list<size_t>::iterator iter = removePoints.begin(); iter != removePoints.end(); iter++) {
            size_t tmpIndex = *iter;
            size_t curPos = tmpIndex - delCounter;

            // GridIndex
            curPoint = list[curPos];

            // erase point
            delCounter++;
            map.erase(curPoint);
            list.erase(list.begin() + curPos);
          }

          // reset all entries in hash map and build list of remaining
          for (size_t i = 0; i < list.size(); i++) {
            curPoint = list[i];
            remainingPoints.push_back(map[curPoint]);
            map[curPoint] = i;
          }

          // reset the whole grid's leaf property in order
          // to guarantee a consistent grid
          recalcLeafProperty();

          // return indices of "surviver"
          return remainingPoints;
        }

        /**
         * unserializes the grid from a string, algorithmic dimensions are not reseted
         *
         * @param istr the string that contains the data
         */
        void unserialize_noAlgoDims(std::string& istr) {
          std::istringstream istream;
          istream.str(istr);

          parseGridDescription(istream);

          //    for (size_t i = 0; i < DIM; i++)
          //    {
          //      algoDims.push_back(i);
          //    }
        }


        /**
         * serialize the gridstorage into a string
         *
         * @return a string that contains all gridstorage information
         */
        std::string serialize() {
          std::ostringstream ostream;
          this->serialize(ostream);
          return ostream.str();
        }

        /**
         * serialize the gridstorage into a stream
         *
         * @param ostream reference to a stream into that all gridstorage information is written
         */
        void serialize(std::ostream& ostream) {
          DimensionBoundary tempBound;

          // Print version, dimensions and number of gridpoints
          ostream << SERIALIZATION_VERSION << " ";
          ostream << DIM << " ";
          ostream << list.size() << std::endl;

          //If BoundingBox used, write zero
          if (!bUseStretching) {
            ostream << std::scientific << 0 << std::endl;

            // Print the bounding box
            for (size_t i = 0; i < DIM; i++) {
              tempBound = boundingBox->getBoundary(i);
              ostream << std::scientific << tempBound.leftBoundary << " " << tempBound.rightBoundary << " " << tempBound.bDirichletLeft << " " << tempBound.bDirichletRight << " ";
            }

            ostream << std::endl;
          } else {
            //If analytic stretching, print the stretching type
            if (*(stretching->getStretchingMode()) == "analytic") {
              ostream << std::scientific << 1 << std::endl;

              // Print the bounding box
              for (size_t i = 0; i < DIM; i++) {
                tempBound = stretching->getBoundary(i);
                ostream << std::scientific << tempBound.leftBoundary << " " << tempBound.rightBoundary << " " << tempBound.bDirichletLeft << " " << tempBound.bDirichletRight << " ";
              }

              ostream << std::endl;
              Stretching1D str1d;
              int stretchingType = 0;

              /*
               * Write stretching type if
               * id: 1
               * log: 2
               * sinh:3
               */
              for (size_t i = 0; i < DIM; i++) {
                str1d = stretching->getStretching1D(i);

                if (str1d.type == "id") {
                  stretchingType = 1;
                } else if (str1d.type == "log") {
                  stretchingType = 2;
                } else if (str1d.type == "sinh") {
                  stretchingType = 3;
                }

                ostream << std::scientific << stretchingType << " " << str1d.x_0 << " " << str1d.xsi << std::endl;

              }
            }
            //If discrete stretching, print the grid vector
            else if (*(stretching->getStretchingMode()) == "discrete") {
              ostream << std::scientific << 2 << std::endl;

              // Print the bounding box
              for (size_t i = 0; i < DIM; i++) {
                tempBound = stretching->getBoundary(i);
                ostream << std::scientific << tempBound.leftBoundary << " " << tempBound.rightBoundary << " " << tempBound.bDirichletLeft << " " << tempBound.bDirichletRight << " ";
              }

              ostream << std::endl;
              std::vector<double>* vec = stretching->getDiscreteVector(true);
              int* vecLevel = stretching->getDiscreteVectorLevel();

              for (size_t i = 0; i < DIM; i++) {
                ostream << std::scientific << vecLevel[i] << std::endl;

                for (size_t j = 0; j < vec[i].size(); j++) {
                  ostream << std::scientific << vec[i][j] << " ";
                }

                ostream << std::endl;
              }
            }
          }

          // print the coordinates of the grid points
          for (grid_list_iterator iter = list.begin(); iter != list.end(); iter++) {
            (*iter)->serialize(ostream);
          }
        }

        /**
        * serialize the gridstorage's gridpoints into a stream
        *
        * @return returns the string that contains all gridpoint information
        */
        std::string toString() {
          std::ostringstream ostream;
          this->toString(ostream);
          return ostream.str();
        }

        /**
         * serialize the gridstorage's gridpoints into a stream
         *
         * @param stream reference to a stream into that all gridpoint information is written
         */
        void toString(std::ostream& stream) {
          stream << "[";
          int i = 0;
          grid_map_iterator iter;

          for (iter = map.begin(); iter != map.end(); iter++, i++) {
            if (i != 0) {
              stream << ",";
            }

            stream << " ";
            iter->first->toString(stream);
            stream << " -> " << iter->second;
          }

          stream << " ]";
        }

        /**
         * gets the size of the hashmap
         *
         * @return returns the size of the hashmap
         */
        size_t size() const {
          return map.size();
        }

        /**
         * gets the number of inner grid points
         *
         * @return the number of inner grid points
         */
        size_t getNumInnerPoints() {
          size_t innerPoints = 0;

          for (size_t p = 0; p < map.size(); p++) {
            if (list[p]->isInnerPoint())
              innerPoints++;
          }

          return innerPoints;
        }

        /**
         * gets the dimension of the grid
         *
         * @return the dimension of the grid stored in this GridStorage object
         */
        size_t dim() const {
          return DIM;
        }

        /**
         * gets the sequence number for given gridpoint by its index
         *
         * @param index gridindex object
         *
         * @return sequence number of the given index
         */
        size_t& operator[](index_pointer index) {
          return map[index];
        }

        /**
         * gets the index number for given gridpoint by its sequence number
         *
         * @param seq the sequence number of the index
         *
         * @return gridindex object (reference)
         */
        index_pointer& operator[](size_t seq) {
          return list[seq];
        }

        /**
         * gets the index number for given gridpoint by its sequence number
         *
         * @param seq the sequence number of the index
         *
         * @return gridindex object (pointer)
         */
        GIT* get(size_t seq) {
          return list[seq];
        }

        /**
         * insert a new index into map
         *
         * @param index reference to the index that should be inserted
         *
         * @return
         */
        size_t insert(index_type& index) {
          index_pointer insert = new GIT(&index);
          list.push_back(insert);
          return (map[insert] = this->seq() - 1);
        }

        /**
         * updates an already stored index
         *
         * @param index reference to the index that should be updated
         * @param pos position where the index should be stored
         */
        void update(index_type& index, size_t pos) {
          if (pos < seq()) {
            // Remove old element at pos
            index_pointer del = list[pos];
            map.erase(del);
            delete del;
            // Insert update
            index_pointer insert = new GIT(&index);
            list[pos] = insert;
            map[insert] = pos;
          }
        }


        /**
         * This methods removes the gridpoint added last. Use with coution, only needed for
         * expanding the grid because of the shadow-storage of prewavelets. Please refer to the
         * Prewavelet grid for further description of the shadow storage.
         *
         */
        void deleteLast() {
          index_pointer del = list.back();
          map.erase(del);
          list.pop_back();
          delete del;
        }

        /**
         * creates a pointer to index from a reference to index by creating
         * a new instance of a index object
         *
         * @param index address of index object
         *
         * @return pointer to new index object
         */
        index_pointer create(index_type& index) {
          index_pointer insert = new GIT(&index);
          return insert;
        }

        /**
         * removes an index from gridstorage
         *
         * @param index pointer to index that should be removed
         */
        void destroy(index_pointer index) {
          delete index;
        }

        /**
         * stores a given index in the hashmap
         *
         * @param index pointer to index that should be stored
         *
         * @return sequence number
         */
        unsigned int store(index_pointer index) {
          list.push_back(index);
          return static_cast<unsigned int>(map[index] = static_cast<unsigned int>(this->seq() - 1));
        }

        /**
         * sets the iterator to a given index
         *
         * @param index the index to which the cursor should be moved
         */
        grid_map_iterator find(index_pointer index) {
          return map.find(index);
        }

        /**
         * set iterator to the first position in the map
         */
        grid_map_iterator begin() {
          return map.begin();
        }

        /**
         * sets the iterator to last position in the map
         */
        grid_map_iterator end() {
          return map.end();
        }

        /**
         * Tests if index is in the storage
         *
         * @param index pointer to index that should be tested
         *
         * @return true if the index is in the storage
         */
        bool has_key(GIT* index) {
          return map.find(index) != map.end();
        }

        /**
         * Sets the index to the left level zero parent
         *
         * @param index pointer to index the should be modified
         * @param dim the dimension in which the modification is taken place
         */
        void left_levelzero(GIT* index, size_t dim) {
          typename index_type::level_type l;
          typename index_type::index_type i;
          index->get(dim, l, i);
          index->set(dim, 0, 0);
        }

        /**
         * Sets the index to the right level zero parent
         *
         * @param index pointer to index the should be modified
         * @param dim the dimension in which the modification is taken place
         */
        void right_levelzero(GIT* index, size_t dim) {
          typename index_type::level_type l;
          typename index_type::index_type i;
          index->get(dim, l, i);
          index->set(dim, 0, 1);
        }

        /**
         * Sets the index to the left child
         *
         * @param index pointer to index the should be modified
         * @param dim the dimension in which the modification is taken place
         */
        void left_child(GIT* index, size_t dim) {
          typename index_type::level_type l;
          typename index_type::index_type i;
          index->get(dim, l, i);
          index->set(dim, l + 1, 2 * i - 1);
        }

        /**
         * Sets the index to the right child
         *
         * @param index pointer to index the should be modified
         * @param dim the dimension in which the modification is taken place
         */
        void right_child(GIT* index, size_t dim) {
          typename index_type::level_type l;
          typename index_type::index_type i;
          index->get(dim, l, i);
          index->set(dim, l + 1, 2 * i + 1);
        }

        /**
         * Resets the index to the top level in direction d
         *
         * @param index pointer to index the should be modified
         * @param d the dimension in which the modification is taken place
         */
        void top(GIT* index, size_t d) {
          index->set(d, 1, 1);
        }

        /**
         * Gets the seq number for index
         *
         * @param index pointer to index which sequence number should be determined
         *
         * @return the seq number for index
         */
        size_t seq(GIT* index) {
          grid_map_iterator iter = map.find(index);

          if (iter != map.end()) {
            return iter->second;
          } else {
            return map.size() + 1;
          }
        }

        /**
         * Tests if seq number does not point to a valid grid index
         *
         * @param s sequence number that should be tested
         *
         * @return true if we are not EOF
         */
        bool end(size_t s) {
          return s > map.size();
        }

        /**
         * returns the algorithmic dimensions (the dimensions in which the Up Down
         * operations should be applied)
         *
         * @return the algorithmic dimensions
         */
        std::vector<size_t> getAlgorithmicDimensions() {
          return algoDims;
        }

        /**
         * sets the algorithmic dimensions (the dimensions in which the Up Down
         * operations should be applied)
         *
         * @param newAlgoDims std::vector containing the algorithmic dimensions
         */
        void setAlgorithmicDimensions(std::vector<size_t> newAlgoDims) {
          algoDims.clear();

          // throw an exception if there is
          if (newAlgoDims.size() > DIM) {
            throw generation_exception("There are more algorithmic dimensions than real dimensions!");
          }

          for (size_t i = 0; i < newAlgoDims.size(); i++) {
            algoDims.push_back(newAlgoDims[i]);
          }
        }

        /**
         * Recalculates the leaf-property of every grid point.
         * This might be useful in case of a grid unserialization
         */
        void recalcLeafProperty() {
          index_pointer point;
          grid_map_iterator iter;
          size_t current_dim;
          typename index_type::level_type l;
          typename index_type::level_type i;
          bool isLeaf = true;

          // iterate through the grid
          for (iter = map.begin(); iter != map.end(); iter++) {
            point = iter->first;
            isLeaf = true;

            // iterate through the dimensions
            for (current_dim = 0; current_dim < DIM; current_dim++) {
              point->get(current_dim, l, i);

              if (l > 0) {
                // Test left child
                left_child(point, current_dim);
                isLeaf = isLeaf && !has_key(point);

                // restore value for dimension
                point->set(current_dim, l, i);

                // Test right child
                right_child(point, current_dim);
                isLeaf = isLeaf && !has_key(point);
              } else {
                // Test level 0
                point->set(current_dim, 1, 1);
                isLeaf = isLeaf && !has_key(point);
              }

              // restore value for dimension
              point->set(current_dim, l, i);
            }

            point->setLeaf(isLeaf);
          }
        }

        /**
         * get the bounding box of the current grid
         *
         * @return returns a pointer to GridStorage's bounding box
         */
        BoundingBox* getBoundingBox() {
          return boundingBox;
        }

        /**
         * get the stretching bounding box of the current grid
         *
         * @return returns a pointer to GridStorage's bounding box
         */
        Stretching* getStretching() {
          return stretching;
        }

        /**
         * sets the bounding box of the current grid
         *
         * @param bb bounding box to which the GridStorage's pointer is set
         */
        void setBoundingBox(BoundingBox& bb) {
          if (!bUseStretching) {
            delete boundingBox;
          }

          if (bUseStretching) {
            delete stretching;
          }

          boundingBox = new BoundingBox(bb);
          bUseStretching = false;
        }

        /**
         * sets the stretching bounding box of the current grid
         *
         * @param bb stretching to which the GridStorage's pointer is set
         */
        void setStretching(Stretching& bb) {
          if (bUseStretching) {
            delete stretching;
          }

          if (!bUseStretching) {
            delete boundingBox;
          }

          bUseStretching = true;
          stretching = new Stretching(bb);
        }

        /**
         * Converts this storage from AOS (array of structures) to SOA (structure of array)
         * with modification to speed up iterative function evaluation. The Level
         * array won't contain the levels, it contains the level to the power of two
         *
         * @param level DataMatrix to store the grid's level to the power of two
         * @param index DataMatrix to store the grid's indices
         */
        void getLevelIndexArraysForEval(DataMatrix& level, DataMatrix& index) {
          typename index_type::level_type curLevel;
          typename index_type::level_type curIndex;

          // Parallelization may lead to segfaults.... comment on your own risk
          //    #pragma omp parallel
          //    {
          //      #pragma omp for schedule (static) private(curLevel, curIndex)
          for (size_t i = 0; i < list.size(); i++) {
            for (size_t current_dim = 0; current_dim < DIM; current_dim++) {
              (list[i])->get(current_dim, curLevel, curIndex);
              level.set(i, current_dim, static_cast<double>(1 << curLevel));
              index.set(i, current_dim, static_cast<double>(curIndex));
            }
          }

          //    }
        }

        /**
         * Converts this storage from AOS (array of structures) to SOA (structure of array)
         * with modification to speed up iterative function evaluation. The Level
         * array won't contain the levels, it contains the level to the power of two
        *
        * The generated arrays are made in format optimized for minimizing page faults
         *
         * @param level DataMatrix to store the grid's level to the power of two
         * @param index DataMatrix to store the grid's indices
         * @param vectorizationType Vectorization type
         * @param blocking_length parameter for an additional blocking length to avoid TLB misses
         */
        void getLevelIndexArraysForEvalTLBOptimized(DataMatrix& level, DataMatrix& index, sg::parallel::VectorizationType vectorizationType, size_t blocking_length) {
          typename index_type::level_type curLevel;
          typename index_type::level_type curIndex;

          //pad datasets
          sg::parallel::DMVectorizationPaddingAssistant::padDataset(level, vectorizationType);
          sg::parallel::DMVectorizationPaddingAssistant::padDataset(index, vectorizationType);

          level.setAll(0.0);
          index.setAll(0.0);

          //transpose
          level.transpose();
          index.transpose();

          //make optimized for reducing page faults

          double* level_ptr = level.getPointer();
          double* index_ptr = index.getPointer();

          for (size_t i = 0; i < list.size(); i += blocking_length) {
            for (size_t current_dim = 0; current_dim < DIM; current_dim++) {
              for (size_t t = i; t < i + blocking_length; ++t) {
                if (t < list.size()) {
                  (list[t])->get(current_dim, curLevel, curIndex);
                  *level_ptr = static_cast<double>(1 << curLevel);
                  *index_ptr = static_cast<double>(curIndex);
                }

                ++level_ptr;
                ++index_ptr;
              }
            }
          }
        }

        /**
         * Converts this storage from AOS (array of structures) to SOA (structure of array)
         * with modification to speed up iterative function evaluation. The Level
         * array won't contain the levels, it contains the level to the power of two
         *
         * @param level DataMatrix to store the grid's level to the power of two
         * @param index DataMatrix to store the grid's indices
         */
        void getLevelIndexArraysForEval(DataMatrixSP& level, DataMatrixSP& index) {
          typename index_type::level_type curLevel;
          typename index_type::level_type curIndex;

          // Parallelization may lead to segfaults.... comment on your own risk
          //    #pragma omp parallel
          //    {
          //      #pragma omp for schedule (static) private(curLevel, curIndex)
          for (size_t i = 0; i < list.size(); i++) {
            for (size_t current_dim = 0; current_dim < DIM; current_dim++) {
              (list[i])->get(current_dim, curLevel, curIndex);
              level.set(i, current_dim, static_cast<float>(1 << curLevel));
              index.set(i, current_dim, static_cast<float>(curIndex));
            }
          }

          //    }
        }

        /**
         * Converts this storage from AOS (array of structures) to SOA (structure of array)
         * with modification to speed up iterative Laplace Calculations: the level
         * won't contain the levels, it contains 2 to the neagative power of the level.
         *
         * @param level DataMatrix to store the grid's modified level
         */
        void getLevelForIntegral(DataMatrix& level) {
          typename index_type::level_type curLevel;
          typename index_type::level_type curIndex;

          // Parallelization may lead to segfaults.... comment on your own risk
          //    #pragma omp parallel
          //    {
          //      #pragma omp for schedule (static) private(curLevel, curIndex)
          for (size_t i = 0; i < list.size(); i++) {
            for (size_t current_dim = 0; current_dim < DIM; current_dim++) {
              (list[i])->get(current_dim, curLevel, curIndex);
              level.set(i, current_dim, pow(2.0, static_cast<int>(-curLevel)));
            }
          }

          //    }
        }

        /**
         * Converts this storage from AOS (array of structures) to SOA (structure of array)
         * with modification to speed up iterative Laplace Calculations: the level
         * won't contain the levels, it contains 2 to the neagative power of the level.
         * Additional blocking for better TLB usage is provided.
         *
         * @param level DataMatrix to store the grid's modified level
         * @param vectorizationType Vectorization type
         * @param blocking_length parameter for an additional blocking length to avoid TLB misses
         */
        void getLevelForIntegralTLBOptimized(DataMatrix& level, sg::parallel::VectorizationType vectorizationType, size_t blocking_length) {
          typename index_type::level_type curLevel;
          typename index_type::level_type curIndex;

          //pad datasets
          sg::parallel::DMVectorizationPaddingAssistant::padDataset(level, vectorizationType);

          level.setAll(0.0);

          //transpose
          level.transpose();

          //make optimized for reducing page faults

          double* level_ptr = level.getPointer();

          for (size_t i = 0; i < list.size(); i += blocking_length) {
            for (size_t current_dim = 0; current_dim < DIM; current_dim++) {
              for (size_t t = i; t < i + blocking_length; ++t) {
                if (t < list.size()) {
                  (list[t])->get(current_dim, curLevel, curIndex);
                  *level_ptr = pow(2.0, static_cast<int>(-curLevel));
                }

                ++level_ptr;
              }
            }
          }
        }

        /**
         * returns the max. depth in all dimension of the grid
         */
        size_t getMaxLevel() {
          typename index_type::level_type curLevel;
          typename index_type::level_type curIndex;
          typename index_type::level_type maxLevel;

          maxLevel = 0;

          for (size_t i = 0; i < list.size(); i++) {
            for (size_t current_dim = 0; current_dim < DIM; current_dim++) {
              (list[i])->get(current_dim, curLevel, curIndex);

              if (curLevel > maxLevel) {
                maxLevel = curLevel;
              }
            }
          }

          return static_cast<size_t>(maxLevel);
        }

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
        void getLevelIndexMaskArraysForModEval(DataMatrix& level, DataMatrix& index, DataMatrix& mask, DataMatrix& offset) {
          typename index_type::level_type curLevel;
          typename index_type::level_type curIndex;

          for (size_t i = 0; i < list.size(); i++) {
            for (size_t current_dim = 0; current_dim < DIM; current_dim++) {
              (list[i])->get(current_dim, curLevel, curIndex);

              if (curLevel == 1) {
                level.set(i, current_dim, 0.0);
                index.set(i, current_dim, 0.0);
                uint64_t intmask = 0x0000000000000000;
                mask.set(i, current_dim, *reinterpret_cast<double*>(&intmask));
                offset.set(i, current_dim, 1.0);
              } else if (curIndex == 1) {
                level.set(i, current_dim, (-1.0)*static_cast<double>(1 << curLevel));
                index.set(i, current_dim, 0.0);
                uint64_t intmask = 0x0000000000000000;
                mask.set(i, current_dim, *reinterpret_cast<double*>(&intmask));
                offset.set(i, current_dim, 2.0);
              } else if (curIndex == static_cast<typename index_type::level_type>( ( (1 << curLevel) - 1) ) ) {
                level.set(i, current_dim, static_cast<double>(1 << curLevel));
                index.set(i, current_dim, static_cast<double>(curIndex));
                uint64_t intmask = 0x0000000000000000;
                mask.set(i, current_dim, *reinterpret_cast<double*>(&intmask));
                offset.set(i, current_dim, 1.0);
              } else {
                level.set(i, current_dim, static_cast<double>(1 << curLevel));
                index.set(i, current_dim, static_cast<double>(curIndex));
                uint64_t intmask = 0x8000000000000000;
                mask.set(i, current_dim, *reinterpret_cast<double*>(&intmask));
                offset.set(i, current_dim, 1.0);
              }
            }
          }
        }

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
        void getLevelIndexMaskArraysForModEval(DataMatrixSP& level, DataMatrixSP& index, DataMatrixSP& mask, DataMatrixSP& offset) {
          typename index_type::level_type curLevel;
          typename index_type::level_type curIndex;

          for (size_t i = 0; i < list.size(); i++) {
            for (size_t current_dim = 0; current_dim < DIM; current_dim++) {
              (list[i])->get(current_dim, curLevel, curIndex);

              if (curLevel == 1) {
                level.set(i, current_dim, 0.0f);
                index.set(i, current_dim, 0.0f);
                unsigned int intmask = 0x00000000;
                mask.set(i, current_dim, *reinterpret_cast<float*>(&intmask));
                offset.set(i, current_dim, 1.0f);
              } else if (curIndex == 1) {
                level.set(i, current_dim, (-1.0f)*static_cast<float>(1 << curLevel));
                index.set(i, current_dim, 0.0f);
                unsigned int intmask = 0x00000000;
                mask.set(i, current_dim, *reinterpret_cast<float*>(&intmask));
                offset.set(i, current_dim, 2.0f);
              } else if (curIndex == static_cast<typename index_type::level_type>(((1 << curLevel) - 1))) {
                level.set(i, current_dim, static_cast<float>(1 << curLevel));
                index.set(i, current_dim, static_cast<float>(curIndex));
                unsigned int intmask = 0x00000000;
                mask.set(i, current_dim, *reinterpret_cast<float*>(&intmask));
                offset.set(i, current_dim, 1.0f);
              } else {
                level.set(i, current_dim, static_cast<float>(1 << curLevel));
                index.set(i, current_dim, static_cast<float>(curIndex));
                unsigned int intmask = 0x80000000;
                mask.set(i, current_dim, *reinterpret_cast<float*>(&intmask));
                offset.set(i, current_dim, 1.0f);
              }
            }
          }
        }

      protected:
        /**
         * returns the next sequence numbers
         *
         * @return returns the next sequence numbers
         */
        size_t seq() {
          return list.size();
        }

      private:
        /// the dimension of the grid
        size_t DIM;
        /// the gridpoints
        grid_list list;
        /// the indecies of the grid points
        grid_map map;
        /// the grid's bounding box
        BoundingBox* boundingBox;
        /// the grid's stretching
        Stretching* stretching;
        /// algorithmic dimension, these are used in Up/Downs
        std::vector<size_t> algoDims;
        /// Flag to check if stretching or bb used
        bool bUseStretching;

        /**
         * Parses the gird's information (grid points, dimensions, bounding box) from a string stream
         *
         * @param istream the string stream that contains the information
         */
        void parseGridDescription(std::istream& istream) {
          int version;
          istream >> version;

          istream >> DIM;

          size_t num;
          istream >> num;

          // check whether grid was created with a version that is too new
          if (version > SERIALIZATION_VERSION) {
            if (version != 4) {
              std::ostringstream errstream;
              errstream << "Version of serialized grid (" << version << ") is too new. Max. recognized version is " << SERIALIZATION_VERSION << ".";
              throw generation_exception(errstream.str().c_str());
            }
          }

          //no bounding box, generate a trivial one
          if (version == 1 || version == 2) {
            // create a standard bounding box
            boundingBox = new BoundingBox(DIM);
          }

          // read the bounding box
          else if (version == 3 || version == 4) {
            // create a standard bounding box
            boundingBox = new BoundingBox(DIM);
            stretching = NULL;
            bUseStretching = false;
            DimensionBoundary tempBound;

            // reads the bounding box
            for (size_t i = 0; i < DIM; i++) {
              istream >> tempBound.leftBoundary;
              istream >> tempBound.rightBoundary;
              istream >> tempBound.bDirichletLeft;
              istream >> tempBound.bDirichletRight;

              boundingBox->setBoundary(i, tempBound);
            }
          } else if (version == 5) {
            //      std::cout<<"Version 5 parse starts\n";
            int useStretching;
            DimensionBoundary tempBound;
            istream >> useStretching;

            if (useStretching == 0) {
              //BoundingBox

              // create a standard bounding box
              boundingBox = new BoundingBox(DIM);
              stretching = NULL;
              bUseStretching = false;


              // reads the boundary data
              for (size_t i = 0; i < DIM; i++) {
                istream >> tempBound.leftBoundary;
                istream >> tempBound.rightBoundary;
                istream >> tempBound.bDirichletLeft;
                istream >> tempBound.bDirichletRight;

                boundingBox->setBoundary(i, tempBound);
              }
            } else if (useStretching == 1) {
              //Stretching with analytic mode
              boundingBox = NULL;
              bUseStretching = true;
              Stretching1D* str1ds = new Stretching1D [DIM];
              DimensionBoundary* tempBounds = new DimensionBoundary[DIM];

              // reads the boundary data
              for (size_t i = 0; i < DIM; i++) {
                istream >> tempBounds[i].leftBoundary;
                istream >> tempBounds[i].rightBoundary;
                istream >> tempBounds[i].bDirichletLeft;
                istream >> tempBounds[i].bDirichletRight;
              }

              int stretchingType = 0;

              //Reads the 1D stretching data
              for (size_t i = 0; i < DIM; i++) {
                istream >> stretchingType;

                switch (stretchingType) {
                  case 1:
                    str1ds[i].type.assign("id");
                    break;

                  case 2:
                    str1ds[i].type.assign("log");
                    break;

                  case 3:
                    str1ds[i].type.assign("sinh");
                    break;

                  default:
                    std::cout << "Stretching Type Unknown in parseGridDescription\n";
                    break;
                }

                istream >> str1ds[i].x_0;
                istream >> str1ds[i].xsi;
              }

              stretching = new Stretching(DIM, tempBounds, str1ds );
              delete [] tempBounds;
              delete[] str1ds;
            } else if (useStretching == 2) {
              //Stretching with discrete Mode

              boundingBox = NULL;
              bUseStretching = true;

              // reads the boundary data, won't be used.
              for (size_t i = 0; i < DIM; i++) {
                istream >> tempBound.leftBoundary;
                istream >> tempBound.rightBoundary;
                istream >> tempBound.bDirichletLeft;
                istream >> tempBound.bDirichletRight;
              }

              int discreteLevel = 0;
              int vectorLength = 0;
              std::vector<double>* vec = new std::vector<double>[DIM];

              for (size_t i = 0; i < DIM; i++) {
                istream >> discreteLevel;
                vectorLength = static_cast<int>(pow(2.0, discreteLevel)) + 1;
                vec[i] = std::vector<double>(vectorLength, 0);

                for (int j = 0; j < vectorLength; j++) {
                  istream >> vec[i][j];
                }
              }

              stretching = new Stretching(DIM, vec);
              delete[] vec;
            } else {
              std::cout << "Unknown Container Id Given in parseGridDescription\n";
            }
          }


          for (size_t i = 0; i < num; i++) {
            index_pointer index = new GIT(istream, version);
            list.push_back(index);
            map[index] = i;
          }

          // set's the grid point's leaf information which is not saved in version 1
          if (version == 1 || version == 4) {
            recalcLeafProperty();
          }
        }
    };

  }
}

#endif /* HASHGRIDSTORAGE_HPP */
