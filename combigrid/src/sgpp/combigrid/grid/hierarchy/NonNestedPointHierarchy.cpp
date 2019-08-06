// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/grid/hierarchy/NonNestedPointHierarchy.hpp>
#include <vector>

namespace sgpp {
namespace combigrid {

NonNestedPointHierarchy::NonNestedPointHierarchy(
    std::shared_ptr<AbstractPointDistribution> pointDistribution,
    std::shared_ptr<AbstractPointOrdering> pointOrdering)
    : numPointsPerLevel(),
      points(),
      permutationIterators(),
      pointDistribution(pointDistribution),
      pointOrdering(pointOrdering) {}

NonNestedPointHierarchy::~NonNestedPointHierarchy() {}

double NonNestedPointHierarchy::getPoint(size_t level, size_t index) {
  return computePoints(level)[index];
}

size_t NonNestedPointHierarchy::getNumPoints(size_t level) {
  while (level >= numPointsPerLevel.size()) {
    size_t currentLevel = numPointsPerLevel.size();
    numPointsPerLevel.push_back(pointOrdering->numPoints(currentLevel));
  }

  return numPointsPerLevel[level];
}

bool NonNestedPointHierarchy::isNested() { return false; }

std::vector<double> &NonNestedPointHierarchy::computePoints(size_t level) {
  while (level >= points.size()) {
    points.push_back(std::vector<double>());
  }

  auto &pointLevel = points[level];
  size_t numPoints = getNumPoints(level);

  while (numPoints > pointLevel.size()) {
    size_t convertedIndex = pointOrdering->convertIndex(level, numPoints, pointLevel.size());
    pointLevel.push_back(pointDistribution->compute(numPoints, convertedIndex));
  }

  return pointLevel;
}

std::vector<double> NonNestedPointHierarchy::getPoints(size_t level, bool sorted) {
  auto &points = computePoints(level);  // could be more than just for this level
  if (!sorted) {
    return points;
  } else {
    size_t numPoints = points.size();
    std::vector<double> result(numPoints);
    auto it = getSortedPermutationIterator(level);
    for (size_t i = 0; i < result.size(); ++i) {
      result[i] = points[it->value()];
      it->moveToNext();
    }

    return result;
  }
}

std::shared_ptr<AbstractPermutationIterator> NonNestedPointHierarchy::getSortedPermutationIterator(
    size_t level) {
  while (level >= permutationIterators.size()) {
    size_t currentLevel = permutationIterators.size();
    size_t numPoints = getNumPoints(currentLevel);
    auto &points = computePoints(currentLevel);
    permutationIterators.push_back(
        pointOrdering->getSortedPermutationIterator(currentLevel, points, numPoints));
  }

  return permutationIterators[level]->clone();
}

} /* namespace combigrid */
} /* namespace sgpp*/
