// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/globaldef.hpp>

#include <set>
#include <map>
#include <vector>

namespace sgpp {
namespace datadriven {

// -------------------------------------------------------------------------------------------

/**
 *
 */
struct HashGridPointCompare {
  bool operator()(const std::shared_ptr<base::HashGridPoint>& lhs,
                  const std::shared_ptr<base::HashGridPoint>& rhs) {
    return lhs->getHash() < rhs->getHash();
  }
};

/**
 *
 */
class OperationMakePositiveCandidateSetAlgorithm {
 public:
  /**
   * Constructor
   *
   * @param maxLevel maximum level for candidate set
   */
  explicit OperationMakePositiveCandidateSetAlgorithm(size_t maxLevel);

  /**
   * Desctructor
   */
  virtual ~OperationMakePositiveCandidateSetAlgorithm();

  /**
   * Load the next candidate set that contains grid points with the currently explored levelsum.
   *
   * @param grid current sparse grid that needs to be extended
   * @param alpha corresponding coefficient vector
   * @param levelSum current levelsum to be explored
   * @param candidates vector that contains the candidate set for the current levelsum
   */
  virtual void nextCandidates(base::Grid& grid, base::DataVector& alpha, size_t levelSum,
                              std::vector<std::shared_ptr<base::HashGridPoint>>& candidates) = 0;

  /**
   * @return number of candidates that have been computed per levelsum
   */
  base::DataVector& numCandidatesPerLevel();

  /**
   * @return number of comparisons that have been computed in total
   */
  size_t costsComputingCandidates();

  /**
   * @return number of comparisons that have been computed per iteration
   */
  base::DataVector& costsComputingCandidatesPerIteration();

  /**
   * @return total number of candidates
   */
  virtual size_t numCandidates() = 0;

  /**
   * @return number of candidates that have been computed per iteration
   */
  virtual base::DataVector& numCandidatesPerIteration() = 0;

  /**
   * Set verbosity level
   *
   * @param pverbose verbosity level
   */
  void setVerbose(bool pverbose);

 protected:
  /**
   * Extract grid points with negative coefficient
   *
   * @param alpha coefficient vector
   * @param negativeGridPoints vector that contains the indices of the grid points with negative
   * coefficient
   * @param tol tolerance for positivity
   */
  void findNodesWithNegativeCoefficients(base::DataVector& alpha,
                                         std::vector<size_t>& negativeGridPoints,
                                         double tol = -1e-14);

  /// iteration counter
  size_t iteration;
  /// maximum full grid level for the candidate set
  size_t maxLevel;
  /// candiddate grid points per level
  base::DataVector gridPointsPerLevel;
  /// comparison costs per iteration
  base::DataVector costsPerIteration;
  /// verbosity level
  bool verbose;
};

// -------------------------------------------------------------------------------------------
class OperationMakePositiveFindIntersectionCandidates
    : public OperationMakePositiveCandidateSetAlgorithm {
 public:
  explicit OperationMakePositiveFindIntersectionCandidates(size_t maxLevel);
  ~OperationMakePositiveFindIntersectionCandidates() override;

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

  virtual void initializeCandidates(base::Grid& grid, std::vector<size_t>& negativeGridPoints);

  void computeIntersection(base::HashGridPoint& gpi, base::HashGridPoint& gpj,
                           base::HashGridPoint& gpintersection);

  static bool compareGridPointsByHash(const std::shared_ptr<base::HashGridPoint>& lhs,
                                      const std::shared_ptr<base::HashGridPoint>& rhs);

  std::unordered_map<size_t, std::shared_ptr<std::vector<std::shared_ptr<base::HashGridPoint>>>>
      intersections;
  std::set<std::shared_ptr<base::HashGridPoint>, HashGridPointCompare> currentIntersections;
  std::set<std::shared_ptr<base::HashGridPoint>, HashGridPointCompare> nextIntersections;
  std::unordered_map<size_t, std::shared_ptr<base::HashGridPoint>> candidates;

  // statistics
  base::DataVector numCandidatesIteration;
};

// -------------------------------------------------------------------------------------------
class OperationMakePositiveFindIntersectionCandidatesJoin
    : public OperationMakePositiveFindIntersectionCandidates {
 public:
  explicit OperationMakePositiveFindIntersectionCandidatesJoin(size_t maxLevel);
  ~OperationMakePositiveFindIntersectionCandidatesJoin() override;

 protected:
  void findIntersections(
      base::Grid& grid, base::DataVector& alpha, size_t levelSum,
      std::unordered_map<size_t, std::shared_ptr<base::HashGridPoint>>& res) override;

  void initializeCandidates(base::Grid& grid, std::vector<size_t>& negativeGridPoints) override;
};

// -------------------------------------------------------------------------------------------
class OperationMakePositiveLoadFullGridCandidates
    : public OperationMakePositiveCandidateSetAlgorithm {
 public:
  explicit OperationMakePositiveLoadFullGridCandidates(size_t maxLevel);
  ~OperationMakePositiveLoadFullGridCandidates() override;

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
  explicit OperationMakePositiveHybridFindIntersectionCandidates(size_t maxLevel);
  ~OperationMakePositiveHybridFindIntersectionCandidates() override;

  void findIntersections(
      base::Grid& grid, base::DataVector& alpha, size_t levelSum,
      std::unordered_map<size_t, std::shared_ptr<base::HashGridPoint>>& res) override;

  void nextCandidates(base::Grid& grid, base::DataVector& alpha, size_t levelSum,
                      std::vector<std::shared_ptr<base::HashGridPoint>>& candidates) override;
};

} /* namespace datadriven */
} /* namespace sgpp */
