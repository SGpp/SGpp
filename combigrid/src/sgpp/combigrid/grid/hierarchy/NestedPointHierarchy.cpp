// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/grid/hierarchy/NestedPointHierarchy.hpp>
#include <vector>

namespace sgpp {
namespace combigrid {

NestedPointHierarchy::NestedPointHierarchy(
    std::shared_ptr<AbstractPointDistribution> pointDistribution,
    std::shared_ptr<AbstractPointOrdering> pointOrdering)
    : numPointsPerLevel(),
      points(),
      permutationIterators(),
      pointDistribution(pointDistribution),
      pointOrdering(pointOrdering) {}

NestedPointHierarchy::~NestedPointHierarchy() {}

double NestedPointHierarchy::getPoint(size_t level, size_t index) {
  // This is a bit hacky, but probably required for CombigridTreeStorage, since
  // it does not store the level in the nested case and thus uses level = 0 for all points
  while (getNumPoints(level) <= index) {
    ++level;
  }

  return computePoints(level)[index];
}

size_t NestedPointHierarchy::getNumPoints(size_t level) {
  while (level >= numPointsPerLevel.size()) {
    size_t currentLevel = numPointsPerLevel.size();
    numPointsPerLevel.push_back(pointOrdering->numPoints(currentLevel));
  }

  return numPointsPerLevel[level];
}

bool NestedPointHierarchy::isNested() { return true; }

std::vector<double>& NestedPointHierarchy::computePoints(size_t level) {
  size_t numPoints = getNumPoints(level);

  while (numPoints > points.size()) {
    size_t convertedIndex = pointOrdering->convertIndex(level, numPoints, points.size());
    points.push_back(pointDistribution->compute(numPoints, convertedIndex));
  }

  return points;
}

std::vector<double> NestedPointHierarchy::getPoints(size_t level, bool sorted) {
  auto& points = computePoints(level);  // could be more than just for this level
  size_t numPoints = getNumPoints(level);
  std::vector<double> result(numPoints);
  if (sorted) {
    auto it = getSortedPermutationIterator(level);
    for (size_t i = 0; i < result.size(); ++i) {
      result[i] = points[it->value()];
      it->moveToNext();
    }
  } else {
    for (size_t i = 0; i < result.size(); ++i) {
      result[i] = points[i];
    }
  }
  return result;
}

std::shared_ptr<AbstractPermutationIterator> NestedPointHierarchy::getSortedPermutationIterator(
    size_t level) {
  while (level >= permutationIterators.size()) {
    size_t currentLevel = permutationIterators.size();
    size_t numPoints = getNumPoints(currentLevel);
    permutationIterators.push_back(
        pointOrdering->getSortedPermutationIterator(currentLevel, points, numPoints));
  }

  return permutationIterators[level]->clone();
}

} /* namespace combigrid */
} /* namespace sgpp*/
