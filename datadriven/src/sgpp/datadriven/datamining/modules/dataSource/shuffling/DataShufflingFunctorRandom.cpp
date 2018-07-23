
/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * DataShufflingFunctorSequential.cpp
 *
 *  Created on: Jul 20, 2018
 *      Author: dominik
 */

#include <sgpp/datadriven/datamining/modules/dataSource/shuffling/DataShufflingFunctorRandom.hpp>

namespace sgpp {
namespace datadriven {

DataShufflingFunctorRandom::DataShufflingFunctorRandom(size_t numSamples, size_t seed) :
    numSamples{numSamples} {
  size_t max = numSamples - 1;

  // Calculate the number of bits needed to store max minus 1
  bitShift = 0;
  while (max != 0) {
    max >>= 1;
    bitShift++;
  }
  bitShift >>= 1;

  // Calculate the bit mask for the right half
  bitMask = (1 << bitShift) - 1;
  std::srand(static_cast<unsigned int>(seed));
  for (size_t i = 0; i < 4; i++) {
    keys[i] = std::rand();
  }
}

size_t DataShufflingFunctorRandom::feistel(size_t idx) {
  // Divide into a left and right part
  size_t left = idx >> bitShift;
  size_t right = idx & bitMask;
  for (size_t i = 0; i < 4; i++) {
    size_t rightNew = left ^ ((hashFunction(right) ^ hashFunction(keys[i])) & bitMask);
    left = right;
    right = rightNew;
  }
  return (right << bitShift) | left;
}

size_t DataShufflingFunctorRandom::operator()(size_t idx) {
  do {
    idx = feistel(idx);
  } while (idx >= numSamples);
  return idx;
}
} /* namespace datadriven */
} /* namespace sgpp */




