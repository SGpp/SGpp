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
#include <set>
#include <utility>
#include <vector>

namespace sgpp {
namespace combigrid {

AdaptiveCombinationGridGenerator::AdaptiveCombinationGridGenerator(
    const std::vector<LevelVector>& levelVectors, const std::vector<double>&& QoIValues,
    std::function<double(double, double)> summationFunction,
    std::shared_ptr<RelevanceCalculator> relevanceCalculator,
    std::shared_ptr<PriorityEstimator> priorityEstimator)
    : summationFunction(summationFunction),
      relevanceCalculator(relevanceCalculator),
      priorityEstimator(priorityEstimator) {
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
    // if the level vector was in the input list, set the corresponding value
    // otherwise, set default value
    auto found = std::find(levelVectors.begin(), levelVectors.end(), level);
    if (found != levelVectors.end()) {
      auto index = distance(levelVectors.begin(), found);
      subspacesAndQoI[level] = QoIValues[index];
    } else {
      subspacesAndQoI[level] = std::numeric_limits<double>::quiet_NaN();
    }
    activeSet.push_back(level);
    adaptLevel(level);
  }
}

AdaptiveCombinationGridGenerator AdaptiveCombinationGridGenerator::fromCombinationGrid(
    const CombinationGrid& combinationGrid, const std::vector<double>&& QoIValues,
    std::function<double(double, double)> summationFunction,
    std::shared_ptr<RelevanceCalculator> relevanceCalculator,
    std::shared_ptr<PriorityEstimator> priorityEstimator) {
  std::vector<LevelVector> levels;

  for (const FullGrid& fullGrid : combinationGrid.getFullGrids()) {
    levels.push_back(fullGrid.getLevel());
  }

  return AdaptiveCombinationGridGenerator(
      levels, std::forward<const std::vector<double>>(QoIValues), summationFunction,
      relevanceCalculator, priorityEstimator);
}

CombinationGrid AdaptiveCombinationGridGenerator::getCombinationGrid(
    const HeterogeneousBasis& basis, bool hasBoundary) const {
  auto oldSet = getOldSet();
  auto coefficients = sgpp::combigrid::getStandardCoefficientsFromLevelSet(oldSet);

  std::vector<LevelVector> oldSetNonzero{};
  // copy only if coefficient is nonzero
  for (size_t i = 0; i < oldSet.size(); ++i) {
    if (coefficients[i] != 0) {
      oldSetNonzero.emplace_back(oldSet[i]);
    }
  }

  return CombinationGrid::fromSubspaces(oldSetNonzero, basis, hasBoundary);
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

double AdaptiveCombinationGridGenerator::getCurrentResult() const {
  double result = 0.;
  auto oldSet = getOldSet();
  auto coefficients = sgpp::combigrid::getStandardCoefficientsFromLevelSet(oldSet);
  // multiply coefficient with QoI and reduce
  for (size_t i = 0; i < oldSet.size(); ++i) {
    if (coefficients[i] != 0.) {
      const auto& value = subspacesAndQoI.at(oldSet[i]);
      assert(!std::isnan(value));
      result = summationFunction(result, coefficients[i] * value);
    }
  }
  return result;
}

std::vector<LevelVector> AdaptiveCombinationGridGenerator::getActiveSet() const {
  return std::vector<LevelVector>(activeSet.begin(), activeSet.end());
}

std::vector<LevelVector> AdaptiveCombinationGridGenerator::getLevels() const {
  auto l = std::vector<LevelVector>(getOldSet());
  auto active = getActiveSet();
  l.insert(l.end(), active.begin(), active.end());
  return l;
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
    if (!std::isnan(subspacesAndQoI.at(lowerHypercube[i]))) {
      level_t hammingDistance = 0;

      for (size_t d = 0; d < levelVector.size(); ++d) {
        hammingDistance += levelVector[d] - lowerHypercube[i][d];
      }
      auto contributionByI = static_cast<double>(subspacesAndQoI.at(lowerHypercube[i]) *
                                                 std::pow(-1, hammingDistance));

      neighborStencilSum = summationFunction(neighborStencilSum, contributionByI);
    }
  }
  assert(!std::isnan(neighborStencilSum));

  return subspacesAndQoI.at(levelVector) + neighborStencilSum;
}

std::vector<double> AdaptiveCombinationGridGenerator::getDeltas(
    const std::vector<LevelVector>& levelVectors) const {
  auto deltas = std::vector<double>();
  deltas.reserve(levelVectors.size());
  for (auto& levelVector : levelVectors) {
    auto delta = getDelta(levelVector);
    deltas.push_back(delta);
  }
  return deltas;
}

std::map<LevelVector, double> AdaptiveCombinationGridGenerator::getPriorities() const {
  std::map<LevelVector, double> priority;
  for (const auto& levelVector : activeSet) {
    std::map<LevelVector, double> deltasOfLowerNeighbors;

    for (size_t d = 0; d < levelVector.size(); ++d) {
      auto neighborLevel = levelVector;
      --neighborLevel[d];
      const bool isInOldSet =
          (std::find(oldSet.begin(), oldSet.end(), neighborLevel) != oldSet.end());
      // check if neighbor exists
      if (isInOldSet) {
        auto delta = getDelta(neighborLevel);
        deltasOfLowerNeighbors[neighborLevel] = delta;
      }
    }

    if (std::isnan(getDelta(levelVector))) {
      priority[levelVector] =
          priorityEstimator->estimatePriority(levelVector, deltasOfLowerNeighbors);
    }
  }
  return priority;
}

std::vector<LevelVector> AdaptiveCombinationGridGenerator::getPriorityQueue() const {
  auto priorities = getPriorities();
  auto comparator = [](std::pair<LevelVector, double> elem1, std::pair<LevelVector, double> elem2) {
    return std::abs(elem1.second) > std::abs(elem2.second);
  };
  std::set<std::pair<LevelVector, double>, decltype(comparator)> orderedSet(
      priorities.begin(), priorities.end(), comparator);

  std::vector<LevelVector> queue;
  queue.reserve(priorities.size());
  std::transform(orderedSet.begin(), orderedSet.end(), back_inserter(queue),
                 [](std::pair<LevelVector, double> const& pair) { return pair.first; });

  return queue;
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

bool AdaptiveCombinationGridGenerator::isAdmissible(const LevelVector& level) const {
  for (size_t d = 0; d < minimumLevelVector.size(); ++d) {
    if (level[d] > minimumLevelVector[d]) {
      LevelVector neighborLevel = level;
      neighborLevel[d] -= 1;
      const bool isInOldSet =
          (std::find(oldSet.begin(), oldSet.end(), neighborLevel) != oldSet.end());

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
