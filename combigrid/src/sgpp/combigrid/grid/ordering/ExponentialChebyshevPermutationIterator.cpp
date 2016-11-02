// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/grid/ordering/ExponentialChebyshevPermutationIterator.hpp>
#include <sgpp/combigrid/utils/Utils.hpp>

namespace sgpp {
namespace combigrid {

ExponentialChebyshevPermutationIterator::ExponentialChebyshevPermutationIterator(
    size_t level, size_t numPoints, size_t currentIndex)
    : currentIndex(currentIndex), level(level), numPoints(numPoints) {}

ExponentialChebyshevPermutationIterator::~ExponentialChebyshevPermutationIterator() {}

void ExponentialChebyshevPermutationIterator::reset() { currentIndex = 0; }

size_t ExponentialChebyshevPermutationIterator::value() {
  size_t reducedIndex = currentIndex;
  size_t reducedLevel = level;

  while (reducedIndex % 3 == 1) {
    --reducedLevel;
    reducedIndex /= 3;
  }

  // preceding points from previous levels + k, if the point is the k-th new point on the level
  return sgpp::combigrid::pow(3, reducedLevel) / 3 + 2 * reducedIndex / 3;
}

void ExponentialChebyshevPermutationIterator::moveToNext() { ++currentIndex; }

std::shared_ptr<AbstractPermutationIterator> ExponentialChebyshevPermutationIterator::clone() {
  return std::shared_ptr<AbstractPermutationIterator>(
      new ExponentialChebyshevPermutationIterator(level, numPoints, currentIndex));
}

} /* namespace combigrid */
} /* namespace sgpp*/
