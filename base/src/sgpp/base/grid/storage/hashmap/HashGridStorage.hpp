// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef HASHGRIDSTORAGE_HPP
#define HASHGRIDSTORAGE_HPP

#include <sgpp/base/exception/generation_exception.hpp>

#include <sgpp/base/grid/storage/hashmap/HashGridPoint.hpp>
#include <sgpp/base/grid/storage/hashmap/SerializationVersion.hpp>

#include <sgpp/base/grid/common/BoundingBox.hpp>
#include <sgpp/base/grid/common/Stretching.hpp>

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataMatrixSP.hpp>

#include <sgpp/globaldef.hpp>

#include <stdint.h>

#include <unordered_map>

#include <exception>
#include <list>
#include <memory>
#include <sstream>
#include <string>
#include <typeinfo>
#include <vector>

namespace sgpp {
namespace base {

class HashGridIterator;

/**
 * Generic hash table based storage of grid points.
 */
class HashGridStorage {
 public:
  /// type of grid points
  typedef HashGridPoint point_type;
  /// pointer to index_type
  typedef HashGridPoint* point_pointer;
  /// pointer to constant index_type
  typedef const HashGridPoint* index_const_pointer;
  /// unordered_map of index_pointers
  typedef std::unordered_map<point_pointer, size_t, HashGridPointPointerHashFunctor,
                             HashGridPointPointerEqualityFunctor>
      grid_map;
  /// iterator of grid_map
  typedef grid_map::iterator grid_map_iterator;
  /// const_iterator of grid_map
  typedef grid_map::const_iterator grid_map_const_iterator;

  /// vector of index_pointers
  typedef std::vector<point_pointer> grid_list;
  /// iterator of grid_list
  typedef grid_list::iterator grid_list_iterator;
  /// const iterator of grid_list
  typedef grid_list::const_iterator grid_list_const_iterator;

  /// iterator for grid points
  typedef HashGridIterator grid_iterator;

  /**
   * Constructor
   *
   * initializes the boundingBox with a trivial cube
   *
   * @param dimension the dimension of the sparse grid
   */
  explicit HashGridStorage(size_t dimension);

  /**
   * Constructor
   *
   * initializes the boundingBox with a reference to another boundingbox
   *
   * @param creationBoundingBox reference to bounding box object that describes the grid boundaries
   */
  explicit HashGridStorage(BoundingBox& creationBoundingBox);

  /**
   * Constructor
   *
   * initializes the stretching with a reference to another stretching
   *
   * @param creationStretching reference to stretching object that describes the grid boundaries and
   * the stretching
   */
  explicit HashGridStorage(Stretching& creationStretching);

  /**
   * Constructor that reads the data from a string
   *
   * @param istr the string that contains the data
   */
  explicit HashGridStorage(std::string& istr);

  /**
   * Constructor that reads the data from an input stream
   *
   * @param istream the inputstream that contains the data
   */
  explicit HashGridStorage(std::istream& istream);

  /**
   * Copy Constructor
   */
  explicit HashGridStorage(HashGridStorage& copyFrom);

  /**
   * Assignment operator
   */
  void operator=(const HashGridStorage& other);

  /**
   * Destructor
   */
  ~HashGridStorage();

  /**
   * deletes all grid points in the storage
   */
  void clear();

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
  void unserializeNoAlgoDims(std::string& istr);

  /**
   * serialize the gridstorage into a string
   *
   * @param version the serialization version of the file
   * @return a string that contains all gridstorage information
   */
  std::string serialize(int version = SERIALIZATION_VERSION) const;

  /**
   * serialize the gridstorage into a stream
   *
   * @param ostream reference to a stream into that all gridstorage information is written
   * @param version the serialization version of the file
   */
  void serialize(std::ostream& ostream, int version = SERIALIZATION_VERSION) const;

  /**
   * serialize the gridstorage's gridpoints into a stream
   *
   * @return returns the string that contains all gridpoint information
   */
  std::string toString() const;

  /**
   * serialize the gridstorage's gridpoints into a stream
   *
   * @param stream reference to a stream into that all gridpoint information is written
   */
  void toString(std::ostream& stream) const;

  /**
   * gets the size of the hashmap
   *
   * @return returns the size of the hashmap
   */
  size_t getSize() const;

  /**
   * gets the number of inner grid points
   *
   * @return the number of inner grid points
   */
  size_t getNumberOfInnerPoints() const;

  /**
   * gets the dimension of the grid
   *
   * @return the dimension of the grid stored in this HashGridStorage object
   */
  size_t getDimension() const;

  /**
   * gets the index number for given gridpoint by its sequence number
   *
   * @param seq the sequence number of the index
   *
   * @return gridpoint object (pointer)
   */
  inline HashGridPoint& operator[](size_t seq) { return *list[seq]; }

  /**
   * gets the index number for given gridpoint by its sequence number
   *
   * @param seq the sequence number of the index
   * @return gridpoint object (constant pointer)
   */
  inline const HashGridPoint& operator[](size_t seq) const { return *list[seq]; }

  /**
   * gets the index number for given gridpoint by its sequence number
   *
   * @param seq the sequence number of the index
   *
   * @return gridpoint object (pointer)
   */
  inline HashGridPoint& getPoint(size_t seq) const { return *list[seq]; }

  /**
   * insert a new index into map
   *
   * @param index reference to the index that should be inserted
   *
   * @return
   */
  size_t insert(const point_type& index);

  /**
   * insert a new index into map including all its ancestors. Boundary points are not added
   *
   * @param index reference to the index that should be inserted
   * @param insertedPoints containing the indices of the new points
   */
  void insert(point_type& index, std::vector<size_t>& insertedPoints);

  /**
   * updates an already stored index
   *
   * @param index reference to the index that should be updated
   * @param pos position where the index should be stored
   */
  void update(point_type& index, size_t pos);

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
  point_pointer create(point_type& index);

  /**
   * removes an index from gridstorage
   *
   * @param index pointer to index that should be removed
   */
  void destroy(point_pointer index);

  /**
   * stores a given index in the hashmap
   *
   * @param index pointer to index that should be stored
   *
   * @return sequence number
   */
  unsigned int store(point_pointer index);

  /**
   * sets the iterator to a given index
   *
   * @param index the index to which the cursor should be moved
   * @return iterator pointing to the index
   */
  grid_map_iterator find(point_pointer index);

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
  bool isContaining(HashGridPoint& index) const;

  /**
   * Gets the seq number for index
   *
   * @param index pointer to index which sequence number should be determined
   *
   * @return the seq number for index
   */
  size_t getSequenceNumber(HashGridPoint& index) const;

  /**
   * Tests if seq number does not point to a valid grid point
   *
   * @param s sequence number that should be tested
   *
   * @return true if we are not EOF
   */
  bool isInvalidSequenceNumber(size_t s);

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
   * @param boundingBox bounding box to which the HashGridStorage's pointer is set
   */
  void setBoundingBox(BoundingBox& boundingBox);

  /**
   * sets the stretching bounding box of the current grid
   *
   * @param stretching stretching to which the HashGridStorage's pointer is set
   */
  void setStretching(Stretching& stretching);

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
   * Converts this storage from AOS (array of structures) to SOA (structure of array)
   * to speed up operations on the position of the grid points.
   *
   * @param coordinates DataMatrix to store the coordinates of the grid points
   */
  void getCoordinateArrays(DataMatrix& coordinates);

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
  void getLevelIndexMaskArraysForModEval(DataMatrix& level, DataMatrix& index, DataMatrix& mask,
                                         DataMatrix& offset);

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

  /**
   * Calculates the coordinate of a given grid point in specific dimension.
   * In contrast to HashGridPoint::getStandardCoordinate, this takes the BoundingBox and
   * Stretching into account.
   *
   * @param point grid point
   * @param d     dimension
   * @return      coordinate of the point in dimension d
   */
  inline double getCoordinate(HashGridPoint point, size_t d) const {
    if ((boundingBox == nullptr) && (stretching == nullptr)) {
      return point.getStandardCoordinate(d);
    } else {
      if (bUseStretching) {
        HashGridPoint::level_type level;
        HashGridPoint::index_type index;

        point.get(d, level, index);
        return stretching->getCoordinate(level, index, d);
      } else {
        return boundingBox->getIntervalWidth(d) * point.getStandardCoordinate(d) +
               boundingBox->getIntervalOffset(d);
      }
    }
  }

  /**
   * Calculates corresponding unit hypercube coordinate of a given point in specific dimension,
   * taking into account the BoundingBox and Stretching.
   */
  inline double getUnitCoordinate(HashGridPoint point, size_t d) const {
    double bbox_point = getCoordinate(point, d);
    if ((boundingBox == nullptr) && (stretching == nullptr)) {
      return bbox_point;
    } else {
      if (bUseStretching) {
        return stretching->transformPointToUnitCube(d, bbox_point);
      } else {
        return boundingBox->transformPointToUnitCube(d, bbox_point);
      }
    }
  }

  /**
   * Calculates the coordinates of a given grid point.
   * In contrast to HashGridPoint::getStandardCoordinates, this takes the BoundingBox and
   * Stretching into account.
   *
   * @param       point         grid point
   * @param[out]  coordinates   vector of coordinates
   */
  void getCoordinates(const HashGridPoint& point, DataVector& coordinates) const;

  /**
   * Calculates the coordinates of a given grid point.
   * In contrast to HashGridPoint::getStandardCoordinates, this takes the BoundingBox and
   * Stretching into account.
   *
   * @param   point   grid point
   * @return          vector of coordinates
   */
  DataVector getCoordinates(const HashGridPoint& point) const;

 private:
  /// the dimension of the grid
  size_t dimension;

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

  /// Flag to check if stretching or boundingBox used
  bool bUseStretching;

  /**
   * Parses the gird's information (grid points, dimensions, bounding box) from a string stream
   *
   * @param istream the string stream that contains the information
   */
  void parseGridDescription(std::istream& istream);
};

HashGridStorage::point_pointer inline HashGridStorage::create(point_type& index) {
  point_pointer insert = new HashGridPoint(index);
  return insert;
}

void inline HashGridStorage::destroy(point_pointer index) { delete index; }

unsigned int inline HashGridStorage::store(point_pointer index) {
  list.push_back(index);
  return static_cast<unsigned int>(map[index] = static_cast<unsigned int>(list.size() - 1));
}

HashGridStorage::grid_map_iterator inline HashGridStorage::find(point_pointer index) {
  return map.find(index);
}

HashGridStorage::grid_map_iterator inline HashGridStorage::begin() { return map.begin(); }

HashGridStorage::grid_map_iterator inline HashGridStorage::end() { return map.end(); }

bool inline HashGridStorage::isContaining(HashGridPoint& index) const {
  return map.find(&index) != map.end();
}

size_t inline HashGridStorage::getSequenceNumber(HashGridPoint& index) const {
  grid_map_const_iterator iter = map.find(&index);

  if (iter != map.end()) {
    return iter->second;
  } else {
    return map.size() + 1;
  }
}

bool inline HashGridStorage::isInvalidSequenceNumber(size_t s) { return s > map.size(); }

std::vector<size_t> inline HashGridStorage::getAlgorithmicDimensions() { return algoDims; }

}  // namespace base
}  // namespace sgpp

#endif /* HASHGRIDSTORAGE_HPP */
