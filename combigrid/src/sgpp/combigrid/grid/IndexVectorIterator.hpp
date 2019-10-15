// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>
#include <sgpp/combigrid/LevelIndexTypes.hpp>
#include <sgpp/combigrid/grid/FullGrid.hpp>

#include <iterator>
#include <vector>

namespace sgpp {
namespace combigrid {

class IndexVectorIterator : public std::iterator<std::random_access_iterator_tag, IndexVector,
                                                 size_t, IndexVector*, IndexVector&> {
 public:
  IndexVectorIterator() : dim(0), indexVector(), minIndex(), maxIndex(),
      numberOfIndexVectors(0), sequenceNumber(0) {
  }

  explicit IndexVectorIterator(const FullGrid& grid) : dim(grid.getDimension()),
      indexVector(dim), minIndex(dim), maxIndex(dim),
      numberOfIndexVectors(dim), sequenceNumber(0) {
    grid.getMinIndex(minIndex);
    grid.getMaxIndex(maxIndex);
    grid.getNumberOfIndexVectors(numberOfIndexVectors);
  }

  IndexVectorIterator(const IndexVector& minIndex, const IndexVector& maxIndex) :
      dim(minIndex.size()), indexVector(dim), minIndex(minIndex), maxIndex(maxIndex),
      numberOfIndexVectors(dim), sequenceNumber(0) {
    for (size_t d = 0; d < dim; d++) {
      numberOfIndexVectors[d] = maxIndex[d] - minIndex[d] + 1;
    }
  }

  IndexVectorIterator(const IndexVectorIterator&) = default;
  IndexVectorIterator& operator=(const IndexVectorIterator&) = default;

  size_t getSequenceNumber() const {
    return sequenceNumber;
  }

  void setSequenceNumber(size_t sequenceNumber) {
    this->sequenceNumber = sequenceNumber;
  }

  IndexVector& operator*() {
    return operator[](sequenceNumber);
  }

  IndexVector* operator->() {
    return &operator[](sequenceNumber);
  }

  IndexVector& operator[](size_t rhs) {
    for (size_t d = 0; d < dim; d++) {
      indexVector[d] = minIndex[d] + static_cast<index_t>(rhs % numberOfIndexVectors[d]);
      rhs /= numberOfIndexVectors[d];
    }

    return indexVector;
  }

  IndexVectorIterator& operator++() {
    sequenceNumber++;
    return *this;
  }

  IndexVectorIterator operator++(int) {
    IndexVectorIterator tmp(*this);
    sequenceNumber++;
    return tmp;
  }

  IndexVectorIterator& operator--() {
    sequenceNumber--;
    return *this;
  }

  IndexVectorIterator operator--(int) {
    IndexVectorIterator tmp(*this);
    sequenceNumber--;
    return tmp;
  }

  IndexVectorIterator operator+(size_t rhs) const {
    IndexVectorIterator tmp(*this);
    tmp.sequenceNumber += rhs;
    return tmp;
  }

  friend IndexVectorIterator operator+(size_t lhs, const IndexVectorIterator& rhs) {
    IndexVectorIterator tmp(rhs);
    tmp.sequenceNumber += lhs;
    return tmp;
  }

  IndexVectorIterator operator-(size_t rhs) const {
    IndexVectorIterator tmp(*this);
    tmp.sequenceNumber -= rhs;
    return tmp;
  }

  size_t operator-(const IndexVectorIterator& other) const {
    return sequenceNumber - other.sequenceNumber;
  }

  IndexVectorIterator& operator+=(size_t rhs) {
    sequenceNumber += rhs;
    return *this;
  }

  IndexVectorIterator& operator-=(size_t rhs) {
    sequenceNumber -= rhs;
    return *this;
  }

  bool operator==(const IndexVectorIterator& other) const {
    return (sequenceNumber == other.sequenceNumber);
  }

  bool operator!=(const IndexVectorIterator& other) const {
    return (sequenceNumber != other.sequenceNumber);
  }

  bool operator<(const IndexVectorIterator& other) const {
    return (sequenceNumber < other.sequenceNumber);
  }

  bool operator>(const IndexVectorIterator& other) const {
    return (sequenceNumber > other.sequenceNumber);
  }

  bool operator<=(const IndexVectorIterator& other) const {
    return (sequenceNumber <= other.sequenceNumber);
  }

  bool operator>=(const IndexVectorIterator& other) const {
    return (sequenceNumber >= other.sequenceNumber);
  }

 protected:
  size_t dim;
  IndexVector indexVector;
  IndexVector minIndex;
  IndexVector maxIndex;
  IndexVector numberOfIndexVectors;
  size_t sequenceNumber;
};

}  // namespace combigrid
}  // namespace sgpp
