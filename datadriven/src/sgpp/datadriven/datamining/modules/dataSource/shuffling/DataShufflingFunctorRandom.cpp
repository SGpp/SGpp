// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/modules/dataSource/shuffling/DataShufflingFunctorRandom.hpp>

#include <iostream>
#include <ctime>
#include <random>

namespace sgpp {
namespace datadriven {

DataShufflingFunctorRandom::DataShufflingFunctorRandom(int64_t seed) :
    numSamples{0}, bitShift{0}, bitMask{0}, seed{seed} {
  if (this->seed == -1) {
    srand(static_cast<unsigned int>(time(nullptr)));
    this->seed = std::rand();
  }
  std::cout << "Seed for random shuffling functor: " << (this->seed) << std::endl;
}

void DataShufflingFunctorRandom::reinitialize(size_t newNumSamples) {
  if (newNumSamples != numSamples && newNumSamples != 0) {
    numSamples = newNumSamples;
    std::mt19937 generator(static_cast<std::mt19937::result_type>(seed));
    size_t max = numSamples - 1;

    // Calculate the number of bits needed to store max minus 1
    bitShift = 0;
    while (max != 0) {
      max >>= 1;
      bitShift++;
    }
    bitShift = -(-bitShift / 2);

    // Calculate the bit mask for the right half
    bitMask = (1 << bitShift) - 1;
    for (size_t i = 0; i < 4; i++) {
      keys[i] = generator();
    }
  }
}

DataShufflingFunctor* DataShufflingFunctorRandom::clone() const {
  return new DataShufflingFunctorRandom{*this};
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

size_t DataShufflingFunctorRandom::operator()(size_t idx, size_t numSamples) {
  reinitialize(numSamples);
  do {
    idx = feistel(idx);
  } while (idx >= numSamples);
  return idx;
}
} /* namespace datadriven */
} /* namespace sgpp */




