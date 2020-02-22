// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/adaptive/AdaptiveCombinationGridGenerator.hpp>

#include <sgpp/base/exception/not_implemented_exception.hpp>

#include <sgpp/combigrid/LevelIndexTypes.hpp>
#include <sgpp/combigrid/grid/CombinationGrid.hpp>

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
    std::unique_ptr<RelevanceCalculator> relevanceCalculator,
    std::unique_ptr<PriorityEstimator> priorityEstimator)
    : relevanceCalculator_(std::move(relevanceCalculator)),
      priorityEstimator_(std::move(priorityEstimator)) {
  assert(levelVectors.size() > 0);

  // set the minimum level vector
  auto numDimensions = levelVectors[0].size();
  minimumLevelVector_ = LevelVector(numDimensions, std::numeric_limits<level_t>::max());
  for (const auto& subspace : levelVectors) {
    for (size_t d = 0; d < numDimensions; ++d) {
      minimumLevelVector_[d] = std::min(subspace[d], minimumLevelVector_[d]);
    }
  }

  auto downwardClosedLevelSet = makeDownwardClosed(levelVectors, minimumLevelVector_);

  for (const auto& level : downwardClosedLevelSet) {
    subspacesAndQoI_[level] = std::numeric_limits<double>::quiet_NaN();
    activeSet_.push_back(level);
    adaptLevel(level);
  }
}

AdaptiveCombinationGridGenerator AdaptiveCombinationGridGenerator::fromCombinationGrid(
    const CombinationGrid& combinationGrid,
    std::unique_ptr<RelevanceCalculator> relevanceCalculator,
    std::unique_ptr<PriorityEstimator> priorityEstimator) {
  std::vector<LevelVector> subspaces{};
  subspaces.resize(combinationGrid.getFullGrids().size());
  std::transform(combinationGrid.getFullGrids().begin(), combinationGrid.getFullGrids().end(),
                 subspaces.begin(),
                 [](const FullGrid& fg) -> LevelVector { return fg.getLevel(); });
  return AdaptiveCombinationGridGenerator(subspaces, std::move(relevanceCalculator),
                                          std::move(priorityEstimator));
}

CombinationGrid AdaptiveCombinationGridGenerator::getCombinationGrid(
    const HeterogeneousBasis& basis) const {
  bool hasBoundary = std::find(minimumLevelVector_.begin(), minimumLevelVector_.end(), 0) !=
                     minimumLevelVector_.end();
  return CombinationGrid::fromSubspaces(oldSet_, basis, hasBoundary);
}

bool AdaptiveCombinationGridGenerator::adaptNextLevelVector(bool regular) {
  if (regular) {
    throw sgpp::base::not_implemented_exception("Parameter regular not yet implemented!");
  }
  auto relevance = getRelevanceOfActiveSet();
  auto pr = std::max_element(
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
  if (subspacesAndQoI_.find(levelVector) == subspacesAndQoI_.end()) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  double neighborStencilSum = 0.;
  auto levelVectorMinusOne = levelVector;
  for (size_t d = 0; d < levelVectorMinusOne.size(); ++d) {
    auto& l = levelVectorMinusOne[d];
    if (l > minimumLevelVector_[d]) {
      l -= 1;
    }
  }

  auto lowerHypercube = hyperCubeOfLevelVectors(levelVector, levelVectorMinusOne);
  lowerHypercube.pop_back();
  // TODO(pollinta): simplify Hamming distance calculation
  for (size_t i = 0; i < lowerHypercube.size(); ++i) {
    auto hypercubeElement = subspacesAndQoI_.find(lowerHypercube[i]);
    assert(hypercubeElement != subspacesAndQoI_.end());
    auto hammingDistance = 0;
    for (size_t d = 0; d < levelVector.size(); ++d) {
      hammingDistance += levelVector[d] - lowerHypercube[i][d];
    }
    neighborStencilSum += subspacesAndQoI_.at(lowerHypercube[i]) * std::pow(-1, hammingDistance);
  }
  return subspacesAndQoI_.at(levelVector) - neighborStencilSum;
}

std::map<LevelVector, double> AdaptiveCombinationGridGenerator::getPriorityQueue() const {
  throw sgpp::base::not_implemented_exception("getPriorityQueue() not yet implemented!");
}

std::map<LevelVector, double> AdaptiveCombinationGridGenerator::getRelevanceOfActiveSet() const {
  std::map<LevelVector, double> relevance{};
  for (const auto& levelVector : activeSet_) {
    const auto delta = getDelta(levelVector);
    if (!std::isnan(delta)) {
      relevance[levelVector] = relevanceCalculator_->calculate(levelVector, delta);
    }
  }
  return relevance;
}

std::map<LevelVector, double>
AdaptiveCombinationGridGenerator::getPrioritiesAndRelevanceOfActiveSet() const {
  auto relevance = getRelevanceOfActiveSet();
  auto priorityQueueAndRelevance = getPriorityQueue();
  priorityQueueAndRelevance.insert(relevance.begin(), relevance.end());
  return priorityQueueAndRelevance;
}

bool AdaptiveCombinationGridGenerator::isAdmissible(const LevelVector& level) const {
  for (size_t d = 0; d < minimumLevelVector_.size(); ++d) {
    if (level[d] > minimumLevelVector_[d]) {
      auto neighborLevel = level;
      neighborLevel[d] -= 1;
      bool isInOldSet = std::find(oldSet_.begin(), oldSet_.end(), neighborLevel) != oldSet_.end();
      if (!isInOldSet) {
        return false;
      }
    }
  }
  return true;
}

void AdaptiveCombinationGridGenerator::addNeighborsToActiveSet(const LevelVector& level) {
  for (size_t d = 0; d < level.size(); ++d) {
    auto neighborLevel = level;
    neighborLevel[d] += 1;
    if (isAdmissible(neighborLevel)) {
      activeSet_.push_back(neighborLevel);
    }
  }
}

void AdaptiveCombinationGridGenerator::adaptLevel(const LevelVector& level) {
  assert(std::find(oldSet_.begin(), oldSet_.end(), level) == oldSet_.end());
  assert(std::find(activeSet_.begin(), activeSet_.end(), level) != activeSet_.end());

  oldSet_.push_back(level);
  activeSet_.remove(level);
  addNeighborsToActiveSet(level);
}

}  // namespace combigrid
}  // namespace sgpp
