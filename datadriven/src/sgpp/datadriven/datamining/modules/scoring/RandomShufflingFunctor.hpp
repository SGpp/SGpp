/* Copyright (C) 2008-today The SG++ project
 *
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * RandomShufflingFunctor.hpp
 *
 *  Created on: 31.07.2016
 *      Author: Michael Lettrich
 */
#pragma once

#include <sgpp/datadriven/datamining/modules/scoring/ShufflingFunctor.hpp>

#include <sgpp/datadriven/tools/Dataset.hpp>

#include <algorithm>
#include <vector>

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
  ShufflingFunctor* clone() const override;

  /**
   * Rearange indices of data samples in a dataset based on a random seed.
   * @param indices vector containing the indices of the data points to be distributed. Vector is
   * modified in place.
   */
  void shuffle(const Dataset& dataset, std::vector<size_t>& indices) override;
};

} /* namespace datadriven */
} /* namespace sgpp */
