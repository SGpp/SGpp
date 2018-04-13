/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * ShufflingFunctor.cpp
 *
 *  Created on: 25.01.2017
 *      Author: Michael Lettrich
 */

#include <sgpp/datadriven/datamining/modules/scoring/ShufflingFunctor.hpp>

namespace sgpp {
namespace datadriven {
ShufflingFunctor::ShufflingFunctor() {
  std::random_device rd;
  seed = rd();
  generator = std::mt19937(seed);
}

int64_t ShufflingFunctor::getSeed() const { return seed; }

void ShufflingFunctor::setSeed(int64_t seed) {
  this->seed = seed;
  generator = std::mt19937(seed);
}
} /* namespace datadriven */
} /* namespace sgpp */
