// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>

#include <sgpp/datadriven/datamining/modules/dataSource/shuffling/DataShufflingFunctor.hpp>

#include <random>

namespace sgpp {
namespace datadriven {

/**
 * A simple shuffling functor for data samples that performs a random shuffling of the data.
 * A Feistel network is used to compute the permutation of [0...N-1] without linear memory
 * complexity.
 */
class DataShufflingFunctorRandom : public DataShufflingFunctor {
 public:
  /**
   * Standard constructor
   * @param seed the seed for the random shuffling (seed == -1 corresponds to a random seed)
   */
  explicit DataShufflingFunctorRandom(int64_t seed = 1337);

  /**
   * Clone pattern.
   * @return identical copy of this instance
   */
  DataShufflingFunctor* clone() const override;

  /**
   * Overload the function-call operator that maps indexes to indexes via a permutation
   * of the entire index set. The permutation used is the identity.
   * @param idx the original index
   * @param numSamples the number of indexes to permute in total
   * @return idx the index after the permutation (simply the input)
   */
  size_t operator()(size_t idx, size_t numSamples) override;

 private:
  /**
   * Next iteration of the feistel network
   *
   */
  size_t feistel(size_t x);

  /**
   * Reinitializes the state with a new number of samples
   */
  void reinitialize(size_t newNumSamples);

  /**
   * The total number of samples in the data source
   */
  size_t numSamples;

  /**
   * The amount of bits that the right part of the key contains
   */
  size_t bitShift;

  /**
   * A bitmask to retrieve the right part of the key (is precomputed to speed up online behaviour)
   */
  size_t bitMask;

  /**
   * Keys for the permutation algorithm
   */
  size_t keys[4];

  /**
   * The hash function to use for the feistel network
   */
  std::hash<size_t> hashFunction;

  /**
   * Store the seed to be able to reinitialize the shuffling
   */
  int64_t seed;
};

} /* namespace datadriven */
} /* namespace sgpp */

