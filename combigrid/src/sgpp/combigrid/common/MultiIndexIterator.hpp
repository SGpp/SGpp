// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_COMMON_MULTIINDEXITERATOR_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_COMMON_MULTIINDEXITERATOR_HPP_

#include <sgpp/combigrid/definitions.hpp>

namespace sgpp {
namespace combigrid {

class MultiIndexIterator {
  MultiIndex index;
  MultiIndex multiBounds;
  bool valid;

 public:
  /**
   * Precondition: all entries in multiBounds are > 0.
   */
  explicit MultiIndexIterator(MultiIndex const &multiBounds)
      : index(multiBounds.size(), 0), multiBounds(multiBounds), valid(true) {}

  void reset() {
    valid = true;
    for (size_t i = 0; i < index.size(); ++i) {
      index[i] = 0;
    }
  }

  bool isValid() { return valid; }

  MultiIndex &value() { return index; }

  size_t indexAt(size_t d) const { return index[d]; }

  MultiIndex getMultiIndex() const { return index; }

  size_t numDimensions() const { return index.size(); }

  int moveToNext() {
    size_t lastDim = index.size() - 1;
    size_t d = lastDim;
    size_t newValue;

    while (true) {
      newValue = (index[d] + 1);
      if (newValue < multiBounds[d]) {
        break;
      }

      if (d == 0) {
        valid = false;
        return -1;
      }

      index[d] = 0;
      --d;
    }

    index[d] = newValue;

    return static_cast<int>(lastDim - d);
  }
};

} /* namespace combigrid */
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_COMMON_MULTIINDEXITERATOR_HPP_ */
