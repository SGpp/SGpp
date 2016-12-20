/*
 * RandomShufflingFunctor.hpp
 *
 *  Created on: 31.07.2016
 *      Author: Michael Lettrich
 */
#pragma once

#include <sgpp/datadriven/datamining/modules/scoring/ShufflingFunctor.hpp>

#include <algorithm>

namespace sgpp {
namespace datadriven {

/**
 * Generate a randomized permutation for a set of indices based on a pseudo random generator. This
 * way sorted samples can be ordered randomly to allow splitting the dataset into testing and
 * training set.
 *
 * ##Usage example
 *
 * vector index        |0|1|2|3|4|5|6
 * --------------------|-|-|-|-|-|-|-
 * index of the sample |0|1|2|3|4|5|6
 *
 * is permuted to e.g.
 *
 * vector index        |0|1|2|3|4|5|6
 * --------------------|-|-|-|-|-|-|-
 * index of the sample |5|1|4|6|3|2|0
 *
 *  indicates where each sample of dataset is moved to in the permuted dataset.
 *
 */
class RandomShufflingFunctor : public ShufflingFunctor {
 public:
  ShufflingFunctor* clone() const override { return new RandomShufflingFunctor{*this}; }

  /**
   * Rearange indices of data samples in a dataset based on a random seed.
   * @param indices vector containing the indices of the data points to be distributed. Vector is
   * modified in place.
   */
  void shuffle(std::vector<size_t>& indices) override {
    std::shuffle(indices.begin(), indices.end(), generator);
  }
};

} /* namespace datadriven */
} /* namespace sgpp */
