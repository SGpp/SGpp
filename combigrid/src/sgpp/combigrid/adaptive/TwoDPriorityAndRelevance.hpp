// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/combigrid/LevelIndexTypes.hpp>
#include <sgpp/combigrid/adaptive/PriorityEstimator.hpp>
#include <sgpp/combigrid/grid/FullGrid.hpp>

#include <map>
#include <vector>

namespace sgpp {
namespace combigrid {

class RelativeMaxRelevanceCalculator : public RelevanceCalculator<std::vector<double>> {
 public:
  explicit WeightedRelevanceCalculator(
      double weightDeltaInRelationToNumberOfPoints = 0.5,
      FullGrid::LevelOccupancy levelOccupancy = FullGrid::LevelOccupancy::TwoToThePowerOfL);

  // needs to get already the relative delta, that is, each entry divided by the average of the
  // predecessors -- PROBLEM //TODO
  double calculate(const LevelVector& levelVector,
                   std::vector<double> delta, std::vector<double> relativeDelta) const override{
                       return std::max(relativeDelta); //
                   }

 private:
  double weightDeltaInRelationToNumberOfPoints;
  FullGrid::LevelOccupancy levelOccupancy;
};

class RelativeMaxAveragingPriorityEstimator : public PriorityEstimator<std::vector<double>> {
 public:
  AveragingPriorityEstimator(
      FullGrid::LevelOccupancy levelOccupancy = FullGrid::LevelOccupancy::TwoToThePowerOfL);

  double estimatePriority(
      const LevelVector& levelVector,
      const std::map<LevelVector, double>& deltasOfDownwardNeighbors) const override;

 private:
  FullGrid::LevelOccupancy levelOccupancy;
};

}  // namespace combigrid
}  // namespace sgpp
