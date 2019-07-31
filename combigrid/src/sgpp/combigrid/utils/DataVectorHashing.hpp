// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataVector.hpp>

#include <vector>
#include <functional>

namespace sgpp {
namespace combigrid {

/**
 * Helper class used internally as an equality predicate.
 */
class DataVectorEqualTo {
 public:
  bool operator()(base::DataVector const& first, base::DataVector const& second) const {
    if (first.getSize() != second.getSize()) {
      return false;
    }

    for (size_t i = 0; i < first.getSize(); ++i) {
      if (first[i] != second[i]) {
        return false;
      }
    }

    return true;
  }
};

/**
 * Helper class used internally as a hash function for DataVector objects.
 */
class DataVectorHash {
 public:
  size_t operator()(base::DataVector const& vec) const {
    std::hash<double> h;
    size_t result = 0;
    for (size_t i = 0; i < vec.getSize(); ++i) {
      result ^= h(vec[i]) * (i + 1);
    }
    return result;
  }
};

} /* namespace combigrid */
} /* namespace sgpp */
