// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/grid/ordering/ExponentialChebyshevPermutationIterator.hpp>
#include <sgpp/combigrid/grid/ordering/ExponentialChebyshevPointOrdering.hpp>
#include <sgpp/combigrid/utils/Utils.hpp>

#include <vector>

namespace sgpp {
namespace combigrid {

ExponentialChebyshevPointOrdering::~ExponentialChebyshevPointOrdering() {}

size_t ExponentialChebyshevPointOrdering::convertIndex(size_t level, size_t numPoints,
                                                       size_t index) {
  size_t firstOccurrenceLevel = 0;
  size_t cumulativePoints = 1;

  while (index >= cumulativePoints) {
    ++firstOccurrenceLevel;
    cumulativePoints *= 3;
  }

  size_t globalIndexInLevel = (index - cumulativePoints / 3);
  // Some elements of this level do belog to previous levels, so we have to add something to the
  // global index
  size_t localIndexInLevel = globalIndexInLevel + (globalIndexInLevel + 1) / 2;
  // stepwidth of first occurrence level compared to current level
  size_t stepwidth = sgpp::combigrid::pow(3, level - firstOccurrenceLevel);

  // this comes from the geometric sum formula (3^w - 1)/(3-1)
  size_t firstIndex = (stepwidth - 1) / 2;

  return firstIndex + localIndexInLevel * stepwidth;
}

size_t ExponentialChebyshevPointOrdering::numPoints(size_t level) {
  return sgpp::combigrid::pow(3, level);
}

std::shared_ptr<AbstractPermutationIterator>
ExponentialChebyshevPointOrdering::getSortedPermutationIterator(size_t level,
                                                                const std::vector<double>& points,
                                                                size_t numPoints) {
  return std::shared_ptr<AbstractPermutationIterator>(
      new ExponentialChebyshevPermutationIterator(level, numPoints));
}

} /* namespace combigrid */
} /* namespace sgpp*/
