// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/combigrid/LevelIndexTypes.hpp>
#include <sgpp/combigrid/grid/FullGrid.hpp>
#include <sgpp/globaldef.hpp>

#include <functional>
#include <iterator>
#include <numeric>
#include <vector>

namespace sgpp {
namespace combigrid {

/**
 * Iterator over the indices contained in a FullGrid.
 */
class IndexVectorIterator : public std::iterator<std::random_access_iterator_tag, IndexVector,
                                                 size_t, IndexVector*, IndexVector&> {
 public:
  /**
   * Default constructor, corresponds to the zero-dimensional case.
   */
  IndexVectorIterator() : dim(0), indexVector(), minIndex(), maxIndex(),
      numberOfIndexVectors(0), sequenceNumber(0) {
  }

  /**
   * Constructor, sets the iterator to the first grid point in the full grid.
   *
   * @param grid  full grid
   */
  explicit IndexVectorIterator(const FullGrid& grid) : dim(grid.getDimension()),
      indexVector(dim), minIndex(dim), maxIndex(dim),
      numberOfIndexVectors(dim), sequenceNumber(0) {
    grid.getMinIndex(minIndex);
    grid.getMaxIndex(maxIndex);
    grid.getNumberOfIndexVectors(numberOfIndexVectors);
  }

  /**
   * Constructor, sets the iterator to the first grid point in the full grid.
   *
   * @param minIndex  vector of minimum 1D indices
   * @param maxIndex  vector of maximum 1D indices
   */
  IndexVectorIterator(const IndexVector& minIndex, const IndexVector& maxIndex) :
      dim(minIndex.size()), indexVector(dim), minIndex(minIndex), maxIndex(maxIndex),
      numberOfIndexVectors(dim), sequenceNumber(0) {
    for (size_t d = 0; d < dim; d++) {
      numberOfIndexVectors[d] = maxIndex[d] - minIndex[d] + 1;
    }
  }

  IndexVectorIterator(const IndexVectorIterator&) = default;
  IndexVectorIterator& operator=(const IndexVectorIterator&) = default;

  /**
   * @return sequence number of current index
   */
  size_t getSequenceNumber() const {
    return sequenceNumber;
  }

  /**
   * @param sequenceNumber   sequence number of current index
   */
  void setSequenceNumber(size_t sequenceNumber) {
    this->sequenceNumber = sequenceNumber;
  }

  /**
   * @return reference to current index
   */
  IndexVector& operator*() {
    return operator[](sequenceNumber);
  }

  /**
   * @return pointer to current index
   */
  IndexVector* operator->() {
    return &operator[](sequenceNumber);
  }

  /**
   * @param rhs   arbitrary sequence number
   * @return reference to index corresponding to given sequence number
   */
  IndexVector& operator[](size_t rhs) {
    for (size_t d = 0; d < dim; d++) {
      indexVector[d] = minIndex[d] + static_cast<index_t>(rhs % numberOfIndexVectors[d]);
      rhs /= numberOfIndexVectors[d];
    }

    return indexVector;
  }

  /**
   * @return iterator before incrementing (selecting the next index)
   */
  IndexVectorIterator& operator++() {
    sequenceNumber++;
    return *this;
  }

  /**
   * @return iterator after incrementing (selecting the next index)
   */
  IndexVectorIterator operator++(int) {
    IndexVectorIterator tmp(*this);
    sequenceNumber++;
    return tmp;
  }

  /**
   * @return iterator before decrementing (selecting the previous index)
   */
  IndexVectorIterator& operator--() {
    sequenceNumber--;
    return *this;
  }

  /**
   * @return iterator after decrementing (selecting the previous index)
   */
  IndexVectorIterator operator--(int) {
    IndexVectorIterator tmp(*this);
    sequenceNumber--;
    return tmp;
  }

  /**
   * @param rhs   right-hand side
   * @return copy of iterator increased by the right-hand side (increase of sequence number)
   */
  IndexVectorIterator operator+(size_t rhs) const {
    IndexVectorIterator tmp(*this);
    tmp.sequenceNumber += rhs;
    return tmp;
  }

  /**
   * @param lhs   left-hand side
   * @param rhs   right-hand side
   * @return copy of right-hand side iterator increased by the left-hand side
   *         (increase of sequence number)
   */
  friend IndexVectorIterator operator+(size_t lhs, const IndexVectorIterator& rhs) {
    IndexVectorIterator tmp(rhs);
    tmp.sequenceNumber += lhs;
    return tmp;
  }

  /**
   * @param rhs   right-hand side
   * @return copy of iterator deccreased by the right-hand side (decrease of sequence number)
   */
  IndexVectorIterator operator-(size_t rhs) const {
    IndexVectorIterator tmp(*this);
    tmp.sequenceNumber -= rhs;
    return tmp;
  }

  /**
   * @param other   other iterator
   * @return difference between this iterator and other iterator (difference of sequence numbers)
   */
  size_t operator-(const IndexVectorIterator& other) const {
    return sequenceNumber - other.sequenceNumber;
  }

  /**
   * @param rhs   right-hand side
   * @return iterator increased by the right-hand side (increase of sequence number)
   */
  IndexVectorIterator& operator+=(size_t rhs) {
    sequenceNumber += rhs;
    return *this;
  }

  /**
   * @param rhs   right-hand side
   * @return iterator decreased by the right-hand side (decrease of sequence number)
   */
  IndexVectorIterator& operator-=(size_t rhs) {
    sequenceNumber -= rhs;
    return *this;
  }

  /**
   * @param other   other iterator
   * @return whether both instances are equal (with respect to sequence numbers)
   */
  bool operator==(const IndexVectorIterator& other) const {
    return (sequenceNumber == other.sequenceNumber);
  }

  /**
   * @param other   other iterator
   * @return whether both instances are inequal (with respect to sequence numbers)
   */
  bool operator!=(const IndexVectorIterator& other) const {
    return (sequenceNumber != other.sequenceNumber);
  }

  /**
   * @param other   other iterator
   * @return whether this instance is smaller than the other instance
   *         (with respect to sequence numbers)
   */
  bool operator<(const IndexVectorIterator& other) const {
    return (sequenceNumber < other.sequenceNumber);
  }

  /**
   * @param other   other iterator
   * @return whether this instance is greater than the other instance
   *         (with respect to sequence numbers)
   */
  bool operator>(const IndexVectorIterator& other) const {
    return (sequenceNumber > other.sequenceNumber);
  }

  /**
   * @param other   other iterator
   * @return whether this instance is smaller than or equal to the other instance
   *         (with respect to sequence numbers)
   */
  bool operator<=(const IndexVectorIterator& other) const {
    return (sequenceNumber <= other.sequenceNumber);
  }

  /**
   * @param other   other iterator
   * @return whether this instance is larger than or equal to the other instance
   *         (with respect to sequence numbers)
   */
  bool operator>=(const IndexVectorIterator& other) const {
    return (sequenceNumber >= other.sequenceNumber);
  }

  /**
   * determines the coordinate of the current index in a given dimension
   * "Standard" means no bounding box (i.e., the domain spanned by minIndex and maxIndex is the unit
   * hypercube) and no stretching (i.e., the points have the standard locations \f$i \cdot
   * 2^{-\ell}\f$).
   *
   * @param d             the dimension in which the coordinate should be calculated
   * @param hasBoundary   whether the grid has boundary points
   *
   * @return the coordinate in the given dimension
   */
  // TODO(pollinta) is this the most sensible place for this functionality?
  inline double getStandardCoordinate(size_t d, bool hasBoundary) {  // todo make const
    auto hInv = hasBoundary ? (maxIndex[d] - minIndex[d]) : (maxIndex[d] - minIndex[d] + 2);
    return static_cast<double>(hasBoundary ? (operator*()[d]) : (operator*()[d]) + 1) /
           static_cast<double>(hInv);
  }

  /**
   * Sets the entries of DataVector coordinates to the coordinates of the gridpoint
   * "Standard" means no bounding box (i.e., the domain spanned by minIndex and maxIndex is the unit
   * hypercube) and no stretching (i.e., the points have the standard locations \f$i \cdot
   * 2^{-\ell}\f$).
   *
   * @param coordinates   the DataVector that should be overwritten with the coordinates
   * @param hasBoundary   whether the grid has boundary points
   */
  void getStandardCoordinates(sgpp::base::DataVector& coordinates, bool hasBoundary = true) {
    coordinates.resize(dim);
    for (size_t d = 0; d < dim; d++) {
      coordinates.set(d, getStandardCoordinate(d, hasBoundary));
    }
  }

  /// check if iterator is at end of range -- needed for swig interface
  bool isAtEnd() {
    static size_t maximumSequenceNumber = std::accumulate(
        numberOfIndexVectors.begin(), numberOfIndexVectors.end(), 1, std::multiplies<index_t>());
    return (dim > 0) && (sequenceNumber == maximumSequenceNumber);
  }

 protected:
  /// dimensionality
  size_t dim;
  /// temporary index vector, overwritten on read access
  IndexVector indexVector;
  /// vector of minimum 1D indices
  IndexVector minIndex;
  /// vector of maximum 1D indices
  IndexVector maxIndex;
  /// number of indices in 1D for all dimensions
  IndexVector numberOfIndexVectors;
  /// current sequence number
  size_t sequenceNumber;
};

}  // namespace combigrid
}  // namespace sgpp
