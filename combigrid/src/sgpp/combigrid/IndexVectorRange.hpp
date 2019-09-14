// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>
#include <sgpp/combigrid/FullGrid.hpp>
#include <sgpp/combigrid/IndexVectorIterator.hpp>
#include <sgpp/combigrid/LevelIndexTypes.hpp>

namespace sgpp {
namespace combigrid {

class IndexVectorRange {
 public:
  IndexVectorRange() : dim(0), grid(), numberOfIndexVectors(0) {
  }

  explicit IndexVectorRange(const FullGrid& grid) : dim(grid.getDimension()), grid(grid),
      numberOfIndexVectors(grid.getNumberOfIndexVectors()) {
  }

  IndexVectorIterator begin() const {
    return IndexVectorIterator::begin(grid);
  }

  IndexVectorIterator end() const {
    return IndexVectorIterator::end(grid);
  }

  size_t find(const IndexVector& index) const {
    size_t result = 0;
    size_t factor = 1;

    for (size_t d = 0; d < dim; d++) {
      result += factor * (index[d] - grid.getMinIndex(d));
      factor *= grid.getNumberOfIndexVectors(d);
    }

    return result;
  }

  const FullGrid& getGrid() const {
    return grid;
  }

  void setGrid(const FullGrid& grid) {
    dim = grid.getDimension();
    this->grid = grid;
    numberOfIndexVectors = grid.getNumberOfIndexVectors();
  }

  void getPoints(base::DataMatrix& points) const {
    const LevelVector& level = grid.getLevel();
    size_t i = 0;
    points.resize(numberOfIndexVectors, dim);

    for (const IndexVector& index : *this) {
      for (size_t d = 0; d < dim; d++) {
        points(i, d) = static_cast<double>(index[d]) /
            static_cast<double>(static_cast<index_t>(1) << level[d]);
      }

      i++;
    }
  }

 protected:
  size_t dim;
  FullGrid grid;
  size_t numberOfIndexVectors;
};

}  // namespace combigrid
}  // namespace sgpp
