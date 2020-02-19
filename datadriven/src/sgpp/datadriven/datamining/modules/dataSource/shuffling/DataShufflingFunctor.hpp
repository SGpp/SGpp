// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>

#include <random>

namespace sgpp {
namespace datadriven {

/**
 * A class to provide functionality to shuffle (reorder) the data samples before the sample
 * provider accesses it. This is neccessary for cross validation and general shuffling of data.
 */
class DataShufflingFunctor {
 public:
  DataShufflingFunctor();
  DataShufflingFunctor(const DataShufflingFunctor& rhs) = default;
  DataShufflingFunctor(DataShufflingFunctor&& rhs) = default;
  DataShufflingFunctor& operator=(const DataShufflingFunctor& rhs) = default;
  DataShufflingFunctor& operator=(DataShufflingFunctor&& rhs) = default;
  virtual ~DataShufflingFunctor() = default;

  /**
   * Polymorphic clone pattern
   * @return deep copy of this object. New object is owned by caller.
   */
  virtual DataShufflingFunctor* clone() const = 0;

  /**
   * Overload the function-call operator that maps indexes to indexes via a permutation
   * of the entire index set.
   * @param idx the original index
   * @param numSamples the number of indexes to permute in total
   * @return idx the index after the permutation
   */
  virtual size_t operator()(size_t idx, size_t numSamples) = 0;
};

} /* namespace datadriven */
} /* namespace sgpp */

