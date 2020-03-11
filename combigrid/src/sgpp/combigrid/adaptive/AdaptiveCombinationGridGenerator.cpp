// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/adaptive/AdaptiveCombinationGridGenerator.hpp>

#include <sgpp/base/exception/not_implemented_exception.hpp>

#include <sgpp/combigrid/LevelIndexTypes.hpp>
#include <sgpp/combigrid/grid/CombinationGrid.hpp>
#include <sgpp/combigrid/tools/LevelVectorTools.hpp>

#include <algorithm>
#include <cassert>
#include <limits>
#include <list>
#include <map>
#include <memory>
#include <numeric>
#include <utility>
#include <vector>

namespace sgpp {
namespace combigrid {

AdaptiveCombinationGridGenerator::AdaptiveCombinationGridGenerator(
    const std::vector<LevelVector>& levelVectors,
    std::function<double(double, double)> summationFunction,
    std::unique_ptr<RelevanceCalculator> relevanceCalculator,
    std::unique_ptr<PriorityEstimator> priorityEstimator)
    : summationFunction(summationFunction),
      relevanceCalculator(std::move(relevanceCalculator)),
      priorityEstimator(std::move(priorityEstimator)) {
  assert(levelVectors.size() > 0);

  // set the minimum level vector
  const size_t numDimensions = levelVectors[0].size();
  minimumLevelVector = LevelVector(numDimensions, std::numeric_limits<level_t>::max());

  for (const LevelVector& subspace : levelVectors) {
    for (size_t d = 0; d < numDimensions; ++d) {
      minimumLevelVector[d] = std::min(subspace[d], minimumLevelVector[d]);
    }
  }

  const std::vector<LevelVector> downwardClosedLevelSet =
      LevelVectorTools::makeDownwardClosed(minimumLevelVector, levelVectors);

  for (const LevelVector& level : downwardClosedLevelSet) {
    subspacesAndQoI[level] = std::numeric_limits<double>::quiet_NaN();
    activeSet.push_back(level);
    adaptLevel(level);
  }
}

AdaptiveCombinationGridGenerator AdaptiveCombinationGridGenerator::fromCombinationGrid(
    const CombinationGrid& combinationGrid, std::function<double(double, double)> summationFunction,
    std::unique_ptr<RelevanceCalculator> relevanceCalculator,
    std::unique_ptr<PriorityEstimator> priorityEstimator) {
  std::vector<LevelVector> subspaces;

  for (const FullGrid& fullGrid : combinationGrid.getFullGrids()) {
    subspaces.push_back(fullGrid.getLevel());
  }

  return AdaptiveCombinationGridGenerator(
      subspaces, summationFunction, std::move(relevanceCalculator), std::move(priorityEstimator));
}

CombinationGrid AdaptiveCombinationGridGenerator::getCombinationGrid(
    const HeterogeneousBasis& basis) const {
  const bool hasBoundary = std::find(minimumLevelVector.begin(), minimumLevelVector.end(), 0) !=
                           minimumLevelVector.end();
  return CombinationGrid::fromSubspaces(oldSet, basis, hasBoundary);
}

bool AdaptiveCombinationGridGenerator::adaptNextLevelVector(bool regular) {
  if (regular) {
    throw sgpp::base::not_implemented_exception("Parameter regular not yet implemented!");
  }

  const std::map<LevelVector, double> relevance = getRelevanceOfActiveSet();
  const auto pr = std::max_element(
      relevance.begin(), relevance.end(),
      [](const MapPairType& m1, const MapPairType& m2) { return m1.second < m2.second; });

  if (pr != relevance.end()) {
    adaptLevel(pr->first);
    return true;
  } else {
    return false;
  }
}

bool AdaptiveCombinationGridGenerator::adaptAllKnown() {
  bool atLeastOneLevelVectorAdded = false;
  bool levelVectorAdded = false;

  do {
    levelVectorAdded = adaptNextLevelVector();
    atLeastOneLevelVectorAdded |= levelVectorAdded;
  } while (levelVectorAdded);

  return atLeastOneLevelVectorAdded;
}

double AdaptiveCombinationGridGenerator::getDelta(const LevelVector& levelVector) const {
  if (subspacesAndQoI.find(levelVector) == subspacesAndQoI.end()) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  double neighborStencilSum = 0.;
  LevelVector levelVectorMinusOne = levelVector;

  for (size_t d = 0; d < levelVectorMinusOne.size(); ++d) {
    if (levelVectorMinusOne[d] > minimumLevelVector[d]) {
      levelVectorMinusOne[d] -= 1;
    }
  }

  std::vector<LevelVector> lowerHypercube =
      LevelVectorTools::generateHyperCube(levelVectorMinusOne, levelVector);
  lowerHypercube.pop_back();

  // TODO(pollinta): simplify Hamming distance calculation
  for (size_t i = 0; i < lowerHypercube.size(); ++i) {
    const auto hypercubeElement = subspacesAndQoI.find(lowerHypercube[i]);
    assert(hypercubeElement != subspacesAndQoI.end());
    level_t hammingDistance = 0;

    for (size_t d = 0; d < levelVector.size(); ++d) {
      hammingDistance += levelVector[d] - lowerHypercube[i][d];
    }

    neighborStencilSum = summationFunction(
        neighborStencilSum,
        static_cast<double>(subspacesAndQoI.at(lowerHypercube[i]) * std::pow(-1, hammingDistance)));
  }

  return subspacesAndQoI.at(levelVector) - neighborStencilSum;
}

std::map<LevelVector, double> AdaptiveCombinationGridGenerator::getPriorityQueue() const {
  throw sgpp::base::not_implemented_exception("getPriorityQueue() not yet implemented!");
}

std::map<LevelVector, double> AdaptiveCombinationGridGenerator::getRelevanceOfActiveSet() const {
  std::map<LevelVector, double> relevance;

  for (const LevelVector& levelVector : activeSet) {
    const double delta = getDelta(levelVector);

    if (!std::isnan(delta)) {
      relevance[levelVector] = relevanceCalculator->calculate(levelVector, delta);
    }
  }

  return relevance;
}

std::map<LevelVector, double>
AdaptiveCombinationGridGenerator::getPrioritiesAndRelevanceOfActiveSet() const {
  const std::map<LevelVector, double> relevance = getRelevanceOfActiveSet();
  std::map<LevelVector, double> priorityQueueAndRelevance = getPriorityQueue();
  priorityQueueAndRelevance.insert(relevance.begin(), relevance.end());
  return priorityQueueAndRelevance;
}

bool AdaptiveCombinationGridGenerator::isAdmissible(const LevelVector& level) const {
  for (size_t d = 0; d < minimumLevelVector.size(); ++d) {
    if (level[d] > minimumLevelVector[d]) {
      LevelVector neighborLevel = level;
      neighborLevel[d] -= 1;
      const bool isInOldSet = (std::find(oldSet.begin(), oldSet.end(), neighborLevel) !=
                               oldSet.end());

      if (!isInOldSet) {
        return false;
      }
    }
  }

  return true;
}

void AdaptiveCombinationGridGenerator::addNeighborsToActiveSet(const LevelVector& level) {
  for (size_t d = 0; d < level.size(); ++d) {
    LevelVector neighborLevel = level;
    neighborLevel[d] += 1;

    if (isAdmissible(neighborLevel)) {
      activeSet.push_back(neighborLevel);
    }
  }
}

void AdaptiveCombinationGridGenerator::adaptLevel(const LevelVector& level) {
  assert(std::find(oldSet.begin(), oldSet.end(), level) == oldSet.end());
  assert(std::find(activeSet.begin(), activeSet.end(), level) != activeSet.end());

  oldSet.push_back(level);
  activeSet.remove(level);
  addNeighborsToActiveSet(level);
}

}  // namespace combigrid
}  // namespace sgpp
