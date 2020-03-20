// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/exception/not_implemented_exception.hpp>
#include <sgpp/combigrid/LevelIndexTypes.hpp>
#include <sgpp/combigrid/grid/FullGrid.hpp>
#include <sgpp/combigrid/tools/IndexVectorIterator.hpp>
#include <sgpp/globaldef.hpp>

#include <vector>

namespace sgpp {
namespace combigrid {

/**
 * Class for iterating over the indices contained in a FullGrid via ranged-based \c for loops.
 * This is basically equivalent to NumPy's and MATLAB's \c meshgrid function, as it loops over
 * the Cartesian product of the sets of 1D indices.
 *
 * The order is as follows: \f$(0, 0, 0, \dotsc, 0)\f$, \f$(1, 0, 0, \dotsc, 0)\f$, ...,
 * \f$(n_1, 0, 0, \dotsc, 0)\f$, \f$(0, 1, 0, \dotsc, 0)\f$, ..., \f$(1, 1, 0, \dotsc, 0)\f$, ...,
 * \f$(n_1, 1, 0, \dotsc, 0)\f$, ..., \f$(n_1, n_2, 0, \dotsc, 0)\f$, etc.
 * (i.e., like enumerating indices of a matrix in column-major order).
 */
class IndexVectorRange {
 public:
  /**
   * Default constructor, corresponds to the zero-dimensional case.
   */
  IndexVectorRange()
      : dim(0), minIndex(), maxIndex(), numberOfIndexVectors(), totalNumberOfIndexVectors(0) {}

  /**
   * Constructor for FullGrid instances.
   *
   * @param grid  full grid
   */
  explicit IndexVectorRange(const FullGrid& grid) : IndexVectorRange() {
    if (grid.getLevelOccupancy() != FullGrid::LevelOccupancy::TwoToThePowerOfL) {
      throw sgpp::base::not_implemented_exception();
    }
    setGrid(grid);
  }

  /**
   * Constructor for custom-defined minimum and maximum indices.
   *
   * @param minIndex  vector of minimum 1D indices
   * @param maxIndex  vector of maximum 1D indices
   */
  IndexVectorRange(const IndexVector& minIndex, const IndexVector& maxIndex)
      : dim(minIndex.size()),
        minIndex(minIndex),
        maxIndex(maxIndex),
        numberOfIndexVectors(dim),
        totalNumberOfIndexVectors(1) {
    for (size_t d = 0; d < dim; d++) {
      numberOfIndexVectors[d] = maxIndex[d] - minIndex[d] + 1;
      totalNumberOfIndexVectors *= numberOfIndexVectors[d];
    }
  }

  /**
   * Returns the beginning of the range, i.e., the IndexVectorIterator that corresponds to
   * the minimum index.
   *
   * @return iterator that corresponds to the minimum index
   */
  IndexVectorIterator begin() const { return IndexVectorIterator(minIndex, maxIndex); }

  /**
   * Returns the end of the range, i.e., the IndexVectorIterator that corresponds to
   * the maximum index.
   *
   * @return iterator that corresponds to the maximum index
   */
  IndexVectorIterator end() const {
    IndexVectorIterator result(minIndex, maxIndex);
    result.setSequenceNumber(totalNumberOfIndexVectors);
    return result;
  }

  /**
   * Finds the given index in the range and returns it position. Only returns reasonable values
   * if the index is actually in the range.
   *
   * @param index   index to find
   * @return sequence number in the range
   */
  size_t find(const IndexVector& index) const {
    size_t result = 0;
    size_t factor = 1;

    for (size_t d = 0; d < dim; d++) {
      result += factor * (index[d] - minIndex[d]);
      factor *= numberOfIndexVectors[d];
    }

    return result;
  }

  /**
   * Sets minimum and maximum index of the range to that of the given FullGrid.
   *
   * @param grid  FullGrid to use
   */
  void setGrid(const FullGrid& grid) {
    dim = grid.getDimension();
    grid.getMinIndex(minIndex);
    grid.getMaxIndex(maxIndex);
    grid.getNumberOfIndexVectors(numberOfIndexVectors);
    totalNumberOfIndexVectors = grid.getNumberOfIndexVectors();
  }

  /**
   * Converts this range to a vector of indices, which will contain all indices that are within
   * this range.
   *
   * @param[out] indices  vector of indices (contents will be overwritten)
   */
  void getIndices(std::vector<IndexVector>& indices) const { indices.assign(begin(), end()); }

  /**
   * Save all grid points of a FullGrid in a DataMatrix.
   *
   * @param[in] grid      full grid
   * @param[out] points   matrix containing the grid points row-by-row after calling
   */
  static void getPoints(const FullGrid& grid, base::DataMatrix& points) {
    IndexVectorRange range(grid);
    const LevelVector& level = grid.getLevel();
    size_t i = 0;
    points.resize(range.totalNumberOfIndexVectors, range.dim);

    if (grid.getLevelOccupancy() != FullGrid::LevelOccupancy::TwoToThePowerOfL) {
      throw sgpp::base::not_implemented_exception();
    }

    for (const IndexVector& index : range) {
      for (size_t d = 0; d < range.dim; d++) {
        points(i, d) = static_cast<double>(index[d]) /
                       static_cast<double>(static_cast<index_t>(1) << level[d]);
      }

      i++;
    }
  }

 protected:
  /// dimensionality
  size_t dim;
  /// vector of minimum 1D indices
  IndexVector minIndex;
  /// vector of maximum 1D indices
  IndexVector maxIndex;
  /// number of indices in 1D for all dimensions
  IndexVector numberOfIndexVectors;
  /// total number of indices
  size_t totalNumberOfIndexVectors;
};

}  // namespace combigrid
}  // namespace sgpp
