// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/grid/ordering/ExponentialNoBoundaryPermutationIterator.hpp>
#include <sgpp/combigrid/grid/ordering/ExponentialNoBoundaryPointOrdering.hpp>
#include <sgpp/combigrid/utils/Utils.hpp>

#include <vector>

namespace sgpp {
namespace combigrid {

ExponentialNoBoundaryPointOrdering::~ExponentialNoBoundaryPointOrdering() {}

size_t ExponentialNoBoundaryPointOrdering::convertIndex(size_t level, size_t numPoints,
                                                        size_t index) {
  size_t firstOccurrenceLevel = 0;
  size_t cumulativePoints = 1;

  while (index >= cumulativePoints) {
    ++firstOccurrenceLevel;
    cumulativePoints = 2 * cumulativePoints + 1;
  }

  size_t globalIndexInLevel = (index - cumulativePoints / 2);
  // Some elements of this level do belog to previous levels, so we have to add something to the
  // global index
  size_t localIndexInLevel = 2 * globalIndexInLevel;
  // stepwidth of first occurrence level compared to current level
  size_t stepwidth = sgpp::combigrid::pow(2, level - firstOccurrenceLevel);

  size_t firstIndex = stepwidth - 1;

  return firstIndex + localIndexInLevel * stepwidth;
}

size_t ExponentialNoBoundaryPointOrdering::numPoints(size_t level) {
  return sgpp::combigrid::pow(2, level + 1) - 1;
}

std::shared_ptr<AbstractPermutationIterator>
ExponentialNoBoundaryPointOrdering::getSortedPermutationIterator(size_t level,
                                                                 const std::vector<double>& points,
                                                                 size_t numPoints) {
  return std::shared_ptr<AbstractPermutationIterator>(
      new ExponentialNoBoundaryPermutationIterator(level, numPoints));
}

} /* namespace combigrid */
} /* namespace sgpp */
