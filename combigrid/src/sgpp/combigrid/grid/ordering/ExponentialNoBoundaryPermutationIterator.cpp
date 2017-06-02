// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/grid/ordering/ExponentialNoBoundaryPermutationIterator.hpp>
#include <sgpp/combigrid/utils/Utils.hpp>

namespace sgpp {
namespace combigrid {

ExponentialNoBoundaryPermutationIterator::ExponentialNoBoundaryPermutationIterator(
    size_t level, size_t numPoints, size_t currentIndex)
    : currentIndex(currentIndex), level(level), numPoints(numPoints) {}

ExponentialNoBoundaryPermutationIterator::~ExponentialNoBoundaryPermutationIterator() {}

void ExponentialNoBoundaryPermutationIterator::reset() { currentIndex = 0; }

size_t ExponentialNoBoundaryPermutationIterator::value() {
  size_t reducedIndex = currentIndex;
  size_t reducedLevel = level;

  while (reducedIndex % 2 == 1) {
    --reducedLevel;
    reducedIndex /= 2;
  }

  // preceding points from previous levels + k, if the point is the k-th new point on the level
  return sgpp::combigrid::pow(2, reducedLevel) - 1 + reducedIndex / 2;
}

void ExponentialNoBoundaryPermutationIterator::moveToNext() { ++currentIndex; }

std::shared_ptr<AbstractPermutationIterator> ExponentialNoBoundaryPermutationIterator::clone() {
  return std::shared_ptr<AbstractPermutationIterator>(
      new ExponentialNoBoundaryPermutationIterator(level, numPoints, currentIndex));
}

} /* namespace combigrid */
} /* namespace sgpp */
