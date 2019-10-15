// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>
#include <sgpp/combigrid/LevelIndexTypes.hpp>
#include <sgpp/combigrid/grid/FullGrid.hpp>
#include <sgpp/combigrid/grid/IndexVectorIterator.hpp>

#include <vector>

namespace sgpp {
namespace combigrid {

class IndexVectorRange {
 public:
  IndexVectorRange() : dim(0), minIndex(), maxIndex(), numberOfIndexVectors(),
      totalNumberOfIndexVectors(0) {
  }

  explicit IndexVectorRange(const FullGrid& grid) : IndexVectorRange() {
    setGrid(grid);
  }

  IndexVectorRange(const IndexVector& minIndex, const IndexVector& maxIndex) :
      dim(minIndex.size()), minIndex(minIndex), maxIndex(maxIndex), numberOfIndexVectors(dim),
      totalNumberOfIndexVectors(1) {
    for (size_t d = 0; d < dim; d++) {
      numberOfIndexVectors[d] = maxIndex[d] - minIndex[d] + 1;
      totalNumberOfIndexVectors *= numberOfIndexVectors[d];
    }
  }

  IndexVectorIterator begin() const {
    return IndexVectorIterator(minIndex, maxIndex);
  }

  IndexVectorIterator end() const {
    IndexVectorIterator result(minIndex, maxIndex);
    result.setSequenceNumber(totalNumberOfIndexVectors);
    return result;
  }

  size_t find(const IndexVector& index) const {
    size_t result = 0;
    size_t factor = 1;

    for (size_t d = 0; d < dim; d++) {
      result += factor * (index[d] - minIndex[d]);
      factor *= numberOfIndexVectors[d];
    }

    return result;
  }

  void setGrid(const FullGrid& grid) {
    dim = grid.getDimension();
    grid.getMinIndex(minIndex);
    grid.getMaxIndex(maxIndex);
    grid.getNumberOfIndexVectors(numberOfIndexVectors);
    totalNumberOfIndexVectors = grid.getNumberOfIndexVectors();
  }

  void getIndices(std::vector<IndexVector>& indices) const {
    indices.assign(begin(), end());
  }

  static void getPoints(const FullGrid& grid, base::DataMatrix& points) {
    IndexVectorRange range(grid);
    const LevelVector& level = grid.getLevel();
    size_t i = 0;
    points.resize(range.totalNumberOfIndexVectors, range.dim);

    for (const IndexVector& index : range) {
      for (size_t d = 0; d < range.dim; d++) {
        points(i, d) = static_cast<double>(index[d]) /
            static_cast<double>(static_cast<index_t>(1) << level[d]);
      }

      i++;
    }
  }

 protected:
  size_t dim;
  IndexVector minIndex;
  IndexVector maxIndex;
  IndexVector numberOfIndexVectors;
  size_t totalNumberOfIndexVectors;
};

}  // namespace combigrid
}  // namespace sgpp
