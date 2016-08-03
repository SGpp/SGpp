// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/grid/ordering/SortedPermutationIterator.hpp>

#include <algorithm>
#include <vector>

namespace sgpp {
namespace combigrid {

SortedPermutationIterator::SortedPermutationIterator(const std::vector<size_t>& permutation,
                                                     size_t currentIndex)
    : permutation(permutation), currentIndex(currentIndex) {}

SortedPermutationIterator::SortedPermutationIterator(const std::vector<double>& points,
                                                     size_t numPoints)
    : permutation(numPoints, 0), currentIndex(0) {
  for (size_t i = 0; i < numPoints; ++i) {
    permutation[i] = i;
  }

  std::sort(permutation.begin(), permutation.begin() + numPoints,
            [&points](size_t i, size_t j) -> bool { return points[i] < points[j]; });
}

SortedPermutationIterator::~SortedPermutationIterator() {}

void SortedPermutationIterator::reset() { currentIndex = 0; }

size_t SortedPermutationIterator::value() { return permutation[currentIndex]; }

void SortedPermutationIterator::moveToNext() { ++currentIndex; }

std::shared_ptr<AbstractPermutationIterator> SortedPermutationIterator::clone() {
  return std::shared_ptr<AbstractPermutationIterator>(
      new SortedPermutationIterator(permutation, currentIndex));
}

} /* namespace combigrid */
} /* namespace sgpp*/
