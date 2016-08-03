// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/grid/ordering/ExponentialChebyshevPointOrdering.hpp>
#include <sgpp/combigrid/grid/ordering/ExponentialChebyshevPermutationIterator.hpp>
#include <sgpp/combigrid/utils/Utils.hpp>

#include <vector>

namespace sgpp {
namespace combigrid {

ExponentialChebyshevPointOrdering::~ExponentialChebyshevPointOrdering() {}

size_t ExponentialChebyshevPointOrdering::convertIndex(size_t level, size_t numPoints,
                                                       size_t index) {
  size_t lastIndex = numPoints - 1;

  if (level == 0) {
    return 0;
  }

  if (index == 0) {
    return lastIndex / 3;
  }
  if (index == 1) {
    return 0;
  }
  if (index == 2) {
    return lastIndex;
  }

  size_t dividedIndex = index - 1;
  size_t resultingLevel = 1;

  // after having subtracted 1, the indices 3, 4 have to be divided once; the indices 5, 6, 7, 8
  // have to be divided twice, etc.
  while (dividedIndex != 1) {
    dividedIndex /= 3;
    ++resultingLevel;
  }

  //  size_t levelHalfPointDistance = (1L << (level - resultingLevel));
  //  size_t indexInLevel = index - ((1L << (resultingLevel - 1)) + 1);
  //
  //  return levelHalfPointDistance + 2 * levelHalfPointDistance * indexInLevel;

  return 0;
}

size_t ExponentialChebyshevPointOrdering::numPoints(size_t level) { return pow(3, level); }

std::shared_ptr<AbstractPermutationIterator>
ExponentialChebyshevPointOrdering::getSortedPermutationIterator(size_t level,
                                                                const std::vector<double>& points,
                                                                size_t numPoints) {
  return std::shared_ptr<AbstractPermutationIterator>(
      new ExponentialChebyshevPermutationIterator(level, numPoints));
}

} /* namespace combigrid */
} /* namespace sgpp*/
