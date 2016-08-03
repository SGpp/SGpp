// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/grid/ordering/IdentityPointOrdering.hpp>
#include <sgpp/combigrid/grid/ordering/SortedPermutationIterator.hpp>

#include <vector>

namespace sgpp {
namespace combigrid {

IdentityPointOrdering::IdentityPointOrdering(std::shared_ptr<AbstractGrowthStrategy> growthStrategy,
                                             bool isSorted)
    : growthStrategy(growthStrategy), isSorted(isSorted) {}

IdentityPointOrdering::~IdentityPointOrdering() {}

size_t IdentityPointOrdering::convertIndex(size_t level, size_t numPoints, size_t index) {
  return index;
}

size_t IdentityPointOrdering::numPoints(size_t level) { return growthStrategy->numPoints(level); }

std::shared_ptr<AbstractPermutationIterator> IdentityPointOrdering::getSortedPermutationIterator(
    size_t level, const std::vector<double>& points, size_t numPoints) {
  if (isSorted) {
    return std::shared_ptr<AbstractPermutationIterator>(nullptr);
  }

  return std::shared_ptr<AbstractPermutationIterator>(
      new SortedPermutationIterator(points, numPoints));
}

} /* namespace combigrid */
} /* namespace sgpp*/
