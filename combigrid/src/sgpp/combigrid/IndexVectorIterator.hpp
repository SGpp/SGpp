// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>
#include <sgpp/combigrid/FullGrid.hpp>
#include <sgpp/combigrid/LevelIndexTypes.hpp>

#include <iterator>
#include <vector>

namespace sgpp {
namespace combigrid {

class IndexVectorIterator : public std::iterator<std::random_access_iterator_tag, IndexVector,
                                                 size_t, IndexVector*, IndexVector&> {
 public:
  explicit IndexVectorIterator(const FullGrid& grid) : dim(grid.getDimension()),
      indexVector(dim), minIndex(dim), maxIndex(dim),
      numberOfIndexVectors(grid.getNumberOfIndexVectors()), sequenceIndex(0) {
    // use move semantics as soon as DataVector supports it
    grid.getMinIndex(minIndex);
    grid.getMaxIndex(maxIndex);
  }

  static IndexVectorIterator begin(const FullGrid& grid) {
    return IndexVectorIterator(grid);
  }

  static IndexVectorIterator end(const FullGrid& grid) {
    IndexVectorIterator result(grid);
    result.setSequenceIndex(grid.getNumberOfIndexVectors());
    return result;
  }

  size_t getSequenceIndex() const {
    return sequenceIndex;
  }

  void setSequenceIndex(size_t sequenceIndex) {
    this->sequenceIndex = sequenceIndex;
  }

  inline void operator++() {
    sequenceIndex++;
  }

  inline size_t operator-(const IndexVectorIterator& other) const {
    return sequenceIndex - other.sequenceIndex;
  }

  inline bool operator==(const IndexVectorIterator& other) const {
    return (sequenceIndex == other.sequenceIndex);
  }

  inline bool operator!=(const IndexVectorIterator& other) const {
    return (sequenceIndex != other.sequenceIndex);
  }

  inline IndexVector& operator*() {
    size_t remainder = sequenceIndex;

    for (size_t d = 0; d < dim; d++) {
      indexVector[d] = minIndex[d] +
          static_cast<index_t>(remainder % (maxIndex[d] - minIndex[d] + 1));
      remainder /= (maxIndex[d] - minIndex[d] + 1);
    }

    return indexVector;
  }

 protected:
  const size_t dim;
  IndexVector indexVector;
  IndexVector minIndex;
  IndexVector maxIndex;
  const size_t numberOfIndexVectors;
  size_t sequenceIndex;
};

}  // namespace combigrid
}  // namespace sgpp
