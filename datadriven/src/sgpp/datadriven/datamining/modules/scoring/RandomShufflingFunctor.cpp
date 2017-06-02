/* Copyright (C) 2008-today The SG++ project
 *
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * RandomShufflingFunctor.cpp
 *
 *  Created on: 25.01.2017
 *      Author: Michael Lettrich
 */

#include <sgpp/datadriven/datamining/modules/scoring/RandomShufflingFunctor.hpp>

#include <sgpp/datadriven/tools/Dataset.hpp>

#include <vector>

namespace sgpp {
namespace datadriven {

ShufflingFunctor* RandomShufflingFunctor::clone() const {
  return new RandomShufflingFunctor{*this};
}

void RandomShufflingFunctor::shuffle(const Dataset& data, std::vector<size_t>& indices) {
  std::shuffle(indices.begin(), indices.end(), generator);
}
} /* namespace datadriven */
} /* namespace sgpp */
