// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/functors/classification/ClassificationRefinementFunctor.hpp>

#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/operation/hash/OperationEval.hpp>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <list>
#include <map>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace sgpp {
namespace datadriven {

ClassificationRefinementFunctor::ClassificationRefinementFunctor(
    std::vector<base::Grid*> grids, std::vector<base::DataVector*> alphas,
    std::vector<double> priors, size_t refinements_num, size_t coarsenings_num,
    bool level_penalize, sgpp::base::AdaptivityThresholdType thresholdType,
    double refinementThreshold, double coarseningThreshold,
    bool coarsenInitialPoints, size_t minimumCoarseningIndex)
    : grids(grids),
      alphas(alphas),
      priors(priors),
      refinements_num(refinements_num),
      coarsenings_num(coarsenings_num),
      level_penalize(level_penalize),
      thresholdType(thresholdType),
      refinementThreshold(refinementThreshold),
      coarseningThreshold(coarseningThreshold),
      coarsenInitialPoints(coarsenInitialPoints),
      total_grid(grids.at(0)->getDimension()),
      minimumCoarseningIndex(minimumCoarseningIndex) {}

double ClassificationRefinementFunctor::operator()(base::GridStorage& storage,
                                                   size_t seq) const {
  return 0.0;
}

std::vector<std::vector<size_t>>
ClassificationRefinementFunctor::adaptAllGrids() {
  // setup total grid
  size_t dim = grids.at(0)->getDimension();
  size_t numClasses = grids.size();
  for (size_t y = 0; y < numClasses; y++) {
    for (size_t seq = 0; seq < grids.at(y)->getSize(); seq++) {
      base::HashGridPoint& gp = grids.at(y)->getStorage().getPoint(seq);
      if (!total_grid.isContaining(gp)) {
        total_grid.insert(gp);
      }
    }
  }

  // initialize neighbors
  base::HashGridPoint root(dim);
  for (size_t j = 0; j < dim; j++) {
    root.getRoot(j);
  }
  root.rehash();
  std::vector<std::pair<base::HashGridPoint, base::HashGridPoint>>
      initialNeighbors(0);
  for (size_t j = 0; j < dim; j++) {
    base::HashGridPoint left(root);
    base::HashGridPoint right(root);
    left.set(j, 0, 0);
    right.set(j, 0, 1);
    initialNeighbors.emplace_back(left, right);
  }
  // start recursion to collect all neighbors
  stepDown(dim, 0, root, initialNeighbors);

  total_grid.recalcLeafProperty();

  // collect coordinates of all grid points
  base::DataMatrix allGridPoints(total_grid.getSize(), dim);
  base::DataVector p(dim);
  for (size_t i = 0; i < total_grid.getSize(); i++) {
    total_grid.getPoint(i).getStandardCoordinates(p);
    allGridPoints.setRow(i, p);
  }

  // evaluate all densities at all grid points
  base::DataMatrix allEvals(total_grid.getSize(), numClasses);
  base::DataVector evals(total_grid.getSize());
  for (size_t y = 0; y < numClasses; y++) {
    sgpp::op_factory::createOperationMultipleEval(*grids.at(y), allGridPoints)
        ->eval(*(alphas.at(y)), evals);
    evals.mult(priors.at(y));
    allEvals.setColumn(y, evals);
  }

  // Obtain dominant classes for each point
  std::vector<size_t> dominantClass(total_grid.getSize());
  for (size_t i = 0; i < total_grid.getSize(); i++) {
    double max = -std::numeric_limits<double>::infinity();
    size_t bestClass = 0;
    for (size_t y = 0; y < numClasses; y++) {
      if (allEvals.get(i, y) > max) {
        max = allEvals.get(i, y);
        bestClass = y;
      }
    }
    dominantClass.at(i) = bestClass;
  }

  // Score neighbor relations
  std::vector<std::multimap<double, std::tuple<size_t, size_t, bool>>>
      classMapsRefinement(numClasses);
  std::vector<std::multimap<double, std::tuple<size_t, size_t, bool>>>
      classMapsCoarsening(numClasses);

  double maxScore = 0.0;

  for (auto const& neighborRel : neighborRels) {
    size_t leafSeqNumber = std::get<0>(neighborRel);
    size_t neighborSeqNumber = std::get<1>(neighborRel);
    size_t classLeaf = dominantClass.at(leafSeqNumber);
    size_t classNeighbor = dominantClass.at(neighborSeqNumber);
    if (classLeaf != classNeighbor) {
      double vLeafClassLeaf = allEvals.get(leafSeqNumber, classLeaf);
      double vLeafClassNeighbor = allEvals.get(leafSeqNumber, classNeighbor);
      double vNeighborClassLeaf = allEvals.get(neighborSeqNumber, classLeaf);
      double vNeighborClassNeighbor =
          allEvals.get(neighborSeqNumber, classNeighbor);
      double score = std::abs((vLeafClassLeaf - vLeafClassNeighbor) -
                              (vNeighborClassLeaf - vNeighborClassNeighbor)) /
                     (1 << std::get<4>(neighborRel));
      std::tuple<size_t, size_t, bool> candidate(
          leafSeqNumber, std::get<2>(neighborRel), std::get<3>(neighborRel));
      classMapsRefinement.at(classLeaf).emplace(score, candidate);
      classMapsRefinement.at(classNeighbor).emplace(score, candidate);

      // check if the grid point of the candidate was already placed in the map
      // (with different
      // neighbor relations)
      insertCoarseningCandidate(classLeaf, classMapsCoarsening, leafSeqNumber,
                                score, candidate);
      insertCoarseningCandidate(classNeighbor, classMapsCoarsening,
                                leafSeqNumber, score, candidate);

      if (score > maxScore) {
        maxScore = score;
      }
    }
  }

  // set the proper threshold based on the max score if a relative threshold is
  // used
  if (thresholdType == sgpp::base::AdaptivityThresholdType::Relative) {
    this->refinementThreshold = this->refinementThreshold * maxScore;
    this->coarseningThreshold = this->coarseningThreshold * maxScore;
  }

  // Insert top refinement candidates into grid
  for (size_t y = 0; y < numClasses; y++) {
    size_t count = 0;
    for (auto it = classMapsRefinement.at(y).rbegin();
         it != classMapsRefinement.at(y).rend(); ++it) {
      if (it->first < this->refinementThreshold) {
        break;
      }
      base::HashGridPoint child(total_grid.getPoint(std::get<0>(it->second)));
      if (std::get<2>(it->second)) {
        child.getLeftChild(std::get<1>(it->second));
      } else {
        child.getRightChild(std::get<1>(it->second));
      }
      std::vector<size_t> insertedPoints;
      grids.at(y)->getStorage().insert(child, insertedPoints);
      count++;
      if (count >= this->refinements_num) {
        break;
      }
    }
  }

  std::vector<std::vector<size_t>> deletedPoints(numClasses);
  // remove top coarsening candidates from the grid
  for (size_t y = 0; y < numClasses; y++) {
    size_t count = 0;
    std::list<size_t> removePoints;

    for (auto it = classMapsCoarsening.at(y).begin();
         it != classMapsCoarsening.at(y).end(); ++it) {
      if (it->first > this->coarseningThreshold) {
        break;
      }

      size_t candidateIndex = std::get<0>(it->second);

      if (!this->coarsenInitialPoints &&
          candidateIndex < this->minimumCoarseningIndex) {
        // do not refine initial grid points
        continue;
      }

      base::HashGridPoint candidate(total_grid.getPoint(candidateIndex));

      removePoints.push_back(candidateIndex);

      count++;
      if (count >= this->coarsenings_num) {
        break;
      }
    }
    std::vector<size_t> idxToDelete{std::begin(removePoints),
                                    std::end(removePoints)};
    deletedPoints.at(y) = idxToDelete;

    if (!removePoints.empty()) {
      grids.at(y)->getStorage().deletePoints(removePoints);
    }
  }

  return deletedPoints;
}

double ClassificationRefinementFunctor::start() const { return 0.0; }

size_t ClassificationRefinementFunctor::getRefinementsNum() const {
  return this->refinements_num;
}

double ClassificationRefinementFunctor::getRefinementThreshold() const {
  return this->refinementThreshold;
}

void ClassificationRefinementFunctor::setGridIndex(size_t grid_index) {
  return;
}

size_t ClassificationRefinementFunctor::getNumGrids() {
  return this->grids.size();
}

void ClassificationRefinementFunctor::preComputeEvaluations() { return; }

void ClassificationRefinementFunctor::stepDown(
    size_t d, size_t minDim, base::HashGridPoint& gp,
    std::vector<std::pair<base::HashGridPoint, base::HashGridPoint>>&
        neighbors) {
  for (size_t j = 0; j < d; j++) {
    base::GridPoint::level_type level;
    base::GridPoint::index_type index;
    gp.get(j, level, index);
    level++;
    index = index << 1;
    base::HashGridPoint left(gp);
    left.set(j, level, index - 1);
    if (!total_grid.isContaining(left)) {
      collectNeighbors(gp, neighbors.at(j).first, j, true);
    } else if (j >= minDim) {
      std::vector<std::pair<base::HashGridPoint, base::HashGridPoint>>
          leftNeighbors(neighbors);
      for (size_t i = 0; i < leftNeighbors.size(); i++) {
        leftNeighbors.at(i).first.set(j, level, index - 1);
        leftNeighbors.at(i).second.set(j, level, index - 1);
      }
      leftNeighbors.at(j) = std::pair<base::HashGridPoint, base::HashGridPoint>(
          neighbors.at(j).first, gp);
      stepDown(d, j, left, leftNeighbors);
    }
    base::HashGridPoint right(gp);
    right.set(j, level, index + 1);
    if (!total_grid.isContaining(right)) {
      collectNeighbors(gp, neighbors.at(j).second, j, false);
    } else if (j >= minDim) {
      std::vector<std::pair<base::HashGridPoint, base::HashGridPoint>>
          rightNeighbors(neighbors);
      for (size_t i = 0; i < rightNeighbors.size(); i++) {
        rightNeighbors.at(i).first.set(j, level, index + 1);
        rightNeighbors.at(i).second.set(j, level, index + 1);
      }
      rightNeighbors.at(j) =
          std::pair<base::HashGridPoint, base::HashGridPoint>(
              gp, neighbors.at(j).second);
      stepDown(d, j, right, rightNeighbors);
    }
  }
}

void ClassificationRefinementFunctor::collectNeighbors(
    base::HashGridPoint leaf, base::HashGridPoint neighbor, size_t dim,
    bool isLeft) {
  size_t leafSeqNumber = total_grid.getSequenceNumber(leaf);
  size_t neighborSeqNumber;
  if (total_grid.isContaining(neighbor)) {
    neighborSeqNumber = total_grid.getSequenceNumber(neighbor);
  } else {
    neighborSeqNumber = total_grid.insert(neighbor);
  }
  neighborRels.emplace_back(leafSeqNumber, neighborSeqNumber, dim, isLeft,
                            leaf.getLevelSum());
}

void ClassificationRefinementFunctor::insertCoarseningCandidate(
    size_t y,
    std::vector<std::multimap<double, std::tuple<size_t, size_t, bool>>>&
        classMapsCoarsening,
    size_t leafSeqNumber, double score,
    std::tuple<size_t, size_t, bool> candidate) {
  for (auto iter = classMapsCoarsening.at(y).begin();
       iter != classMapsCoarsening.at(y).end();) {
    // we are looking for the candidate with the highest score -> only insert
    // new candidate if new
    // score is higher
    if (std::get<0>(iter->second) == leafSeqNumber && score > iter->first) {
      // discard old candidate
      iter = classMapsCoarsening.at(y).erase(iter);
      break;
    } else if (std::get<0>(iter->second) == leafSeqNumber) {
      // discard new candidate, candidate with higher score already exists.
      return;
    }
    ++iter;
  }
  // insert the point if it is not already in the map or the new candidate has a
  // higher score
  classMapsCoarsening.at(y).emplace(score, candidate);
}

}  // namespace datadriven
}  // namespace sgpp
