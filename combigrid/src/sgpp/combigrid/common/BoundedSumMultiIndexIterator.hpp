// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_COMMON_BOUNDEDSUMMULTIINDEXITERATOR_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_COMMON_BOUNDEDSUMMULTIINDEXITERATOR_HPP_

#include <sgpp/combigrid/definitions.hpp>

#include <cstddef>

namespace sgpp {
namespace combigrid {

/**
 * This class can be used to generate multi-indices with bounded 1-norm. For example, using a
 * BoundedSumMultiIndexIterator with dim = 2 and maxIndexSum = 2 will generate the indices
 * (0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (2, 0).
 */
class BoundedSumMultiIndexIterator {
  MultiIndex index;
  size_t maxIndexSum;
  size_t indexSum;
  bool valid;

 public:
  BoundedSumMultiIndexIterator(size_t dim, size_t maxIndexSum)
      : index(dim, 0), maxIndexSum(maxIndexSum), indexSum(0), valid(true) {}

  /**
   * Start again from the beginning of the sequence.
   */
  void reset() {
    valid = true;
    for (size_t i = 0; i < index.size(); ++i) {
      index[i] = 0;
    }
  }

  /**
   * The iterator is valid iff the current value exists (i. e. the end of the sequence has not been
   * reached yet).
   */
  bool isValid() { return valid; }

  /**
   * Returns the current multi-index.
   */
  MultiIndex &value() { return index; }

  /**
   * Returns a particular component (in dimension d, starting from 0) of the current multi-index.
   */
  size_t indexAt(size_t d) const { return index[d]; }

  /**
   * Moves the iterator forward to the next multi-index.
   */
  int moveToNext() {
    size_t lastDim = index.size() - 1;
    size_t d = lastDim;

    while (indexSum >= maxIndexSum) {
      if (d == 0) {
        valid = false;
        return -1;
      }

      indexSum -= index[d];
      index[d] = 0;

      --d;
    }

    ++index[d];
    ++indexSum;

    return static_cast<int>(lastDim - d);
  }
};

} /* namespace combigrid */
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_COMMON_BOUNDEDSUMMULTIINDEXITERATOR_HPP_ */
