// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/globaldef.hpp>

#include <map>
#include <vector>

namespace sgpp {
namespace datadriven {

// -------------------------------------------------------------------------------------------

class OperationMakePositiveCandidateSetAlgorithm {
 public:
  OperationMakePositiveCandidateSetAlgorithm(size_t maxLevel);
  virtual ~OperationMakePositiveCandidateSetAlgorithm();

  virtual void nextCandidates(base::Grid& grid, base::DataVector& alpha, size_t levelSum,
                              std::vector<std::shared_ptr<base::HashGridPoint>>& candidates) = 0;

  // statistics
  base::DataVector& numCandidatesPerLevel();
  size_t costsComputingCandidates();
  base::DataVector& costsComputingCandidatesPerIteration();

  virtual size_t numCandidates() = 0;
  virtual base::DataVector& numCandidatesPerIteration() = 0;

  void setVerbose(bool pverbose);

 protected:
  void findNodesWithNegativeCoefficients(base::DataVector& alpha,
                                         std::vector<size_t>& negativeGridPoints,
                                         double tol = -1e-14);

  size_t iteration;
  size_t maxLevel;
  base::DataVector gridPointsPerLevel;
  base::DataVector costsPerIteration;
  bool verbose;
};

// -------------------------------------------------------------------------------------------
class OperationMakePositiveFindIntersectionCandidates
    : public OperationMakePositiveCandidateSetAlgorithm {
 public:
  OperationMakePositiveFindIntersectionCandidates(size_t maxLevel);
  virtual ~OperationMakePositiveFindIntersectionCandidates();

  void nextCandidates(base::Grid& grid, base::DataVector& alpha, size_t levelSum,
                      std::vector<std::shared_ptr<base::HashGridPoint>>& candidates) override;

  size_t numCandidates() override;
  base::DataVector& numCandidatesPerIteration() override;

 protected:
  bool haveOverlappingSupport(base::HashGridPoint& gpi, base::HashGridPoint& gpj, size_t dim);
  bool haveOverlappingSupport(base::HashGridPoint& gpi, base::HashGridPoint& gpj);

  virtual void findIntersections(
      base::Grid& grid, base::DataVector& alpha, size_t levelSum,
      std::unordered_map<size_t, std::shared_ptr<base::HashGridPoint>>& res);

  void initializeCandidates(base::Grid& grid, std::vector<size_t>& negativeGridPoints);

  void computeIntersection(base::HashGridPoint& gpi, base::HashGridPoint& gpj,
                           base::HashGridPoint& gpintersection);

  static bool compareGridPointsByHash(const std::shared_ptr<base::HashGridPoint>& lhs,
                                      const std::shared_ptr<base::HashGridPoint>& rhs);

  std::unordered_map<size_t, std::shared_ptr<std::vector<std::shared_ptr<base::HashGridPoint>>>>
      intersections;
  std::unordered_map<size_t, std::shared_ptr<base::HashGridPoint>> currentIntersections;
  std::unordered_map<size_t, std::shared_ptr<base::HashGridPoint>> nextIntersections;
  std::unordered_map<size_t, std::shared_ptr<base::HashGridPoint>> candidates;

  // statistics
  base::DataVector numCandidatesIteration;
};

// -------------------------------------------------------------------------------------------
class OperationMakePositiveLoadFullGridCandidates
    : public OperationMakePositiveCandidateSetAlgorithm {
 public:
  OperationMakePositiveLoadFullGridCandidates(size_t maxLevel);
  virtual ~OperationMakePositiveLoadFullGridCandidates();

  void nextCandidates(base::Grid& grid, base::DataVector& alpha, size_t levelSum,
                      std::vector<std::shared_ptr<base::HashGridPoint>>& candidates) override;

  size_t numCandidates() override;
  base::DataVector& numCandidatesPerIteration() override;

 private:
  void initializeFullGrid(base::Grid& grid);

  std::unique_ptr<base::Grid> fullGrid;
  base::DataVector candidatesPerIteration;
};

// -------------------------------------------------------------------------------------------
class OperationMakePositiveHybridFindIntersectionCandidates
    : public OperationMakePositiveFindIntersectionCandidates {
 public:
  OperationMakePositiveHybridFindIntersectionCandidates(size_t maxLevel);
  virtual ~OperationMakePositiveHybridFindIntersectionCandidates();

  void findIntersections(
      base::Grid& grid, base::DataVector& alpha, size_t levelSum,
      std::unordered_map<size_t, std::shared_ptr<base::HashGridPoint>>& res) override;

  void nextCandidates(base::Grid& grid, base::DataVector& alpha, size_t levelSum,
                      std::vector<std::shared_ptr<base::HashGridPoint>>& candidates) override;

 protected:
  size_t overallComparisonCosts;
};

} /* namespace datadriven */
} /* namespace sgpp */
