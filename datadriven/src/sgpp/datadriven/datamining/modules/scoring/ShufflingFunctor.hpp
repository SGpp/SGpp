/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * ShufflingFunctor.hpp
 *
 *  Created on: 31.07.2016
 *      Author: Michael Lettrich
 */

#pragma once

#include <sgpp/globaldef.hpp>

#include <random>
#include <vector>

namespace sgpp {
namespace datadriven {

/**
 * A shuffling functor generates a permutation for a set of indices. Can be used to generate a
 * random or particular ordering of samples from a potentially ordered dataset as training and
 * testing data should be evenly distributed when performing training with testing or
 * cross-validation.
 */
class ShufflingFunctor {
 public:
  ShufflingFunctor();
  ShufflingFunctor(const ShufflingFunctor& rhs) = default;
  ShufflingFunctor(ShufflingFunctor&& rhs) = default;
  ShufflingFunctor& operator=(const ShufflingFunctor& rhs) = default;
  ShufflingFunctor& operator=(ShufflingFunctor&& rhs) = default;
  virtual ~ShufflingFunctor() = default;

  /**
   * Polymorphic clone pattern
   * @return deep copy of this object. New object is owned by caller.
   */
  virtual ShufflingFunctor* clone() const = 0;

  /**
   * Create a permutation from a vector of indices. The indices can then be mapped back to samples.
   * @param indices: vector containing indices. Will permute indices in place.
   */
  virtual void shuffle(std::vector<size_t>& indices) = 0;

  /**
   * Get random seed for randomized operations
   * @return random seed.
   */
  int64_t getSeed() const;

  /**
   * Set random seed for randomized operations
   * @param seed the new random seed.
   */
  void setSeed(int64_t seed);

 protected:
  /**
   * Random seed for a pseudo random number generator. Same seeds should produce identical random
   * numbers. Required as an initialization for the random number generator.
   */
  int64_t seed;

  /**
   * Pseudo random number generator required for shuffling.
   */
  std::mt19937 generator;
};

} /* namespace datadriven */
} /* namespace sgpp */
