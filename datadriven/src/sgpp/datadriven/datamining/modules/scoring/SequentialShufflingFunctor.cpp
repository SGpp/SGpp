/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * SequentialShufflingFunctor.cpp
 *
 *  Created on: 25.01.2017
 *      Author: Michael Lettrich
 */
#include <sgpp/datadriven/datamining/modules/scoring/SequentialShufflingFunctor.hpp>
#include <vector>

namespace sgpp {
namespace datadriven {
ShufflingFunctor* SequentialShufflingFunctor::clone() const {
  return new SequentialShufflingFunctor{*this};
}

void SequentialShufflingFunctor::shuffle(const Dataset& data, std::vector<size_t>& indices) {
  // doesn't do any permutation to provided indices, since we're using the entries sequentially
}
} /* namespace datadriven */
} /* namespace sgpp */
