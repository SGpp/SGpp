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

class OperationMakePositiveCandidateSetAlgorithm {
 public:
  OperationMakePositiveCandidateSetAlgorithm();
  virtual ~OperationMakePositiveCandidateSetAlgorithm();

  virtual void nextCandidates(base::Grid& grid, base::DataVector& alpha, size_t levelSum,
                              std::vector<std::shared_ptr<base::HashGridPoint>>& candidates) = 0;

  virtual size_t numCandidates() = 0;
  void setVerbose(bool pverbose);

 protected:
  void findNodesWithNegativeCoefficients(base::DataVector& alpha,
                                         std::vector<size_t>& negativeGridPoints,
                                         double tol = -1e-14);

  size_t iteration;
  bool verbose;
};
// -------------------------------------------------------------------------------------------
class OperationMakePositiveFindIntersectionCandidates
    : public OperationMakePositiveCandidateSetAlgorithm {
 public:
  OperationMakePositiveFindIntersectionCandidates();
  virtual ~OperationMakePositiveFindIntersectionCandidates();

  void nextCandidates(base::Grid& grid, base::DataVector& alpha, size_t levelSum,
                      std::vector<std::shared_ptr<base::HashGridPoint>>& candidates) override;

  size_t numCandidates() override;

 protected:
  bool haveOverlappingSupport(base::HashGridPoint& gpi, base::HashGridPoint& gpj, size_t dim);
  bool haveOverlappingSupport(base::HashGridPoint& gpi, base::HashGridPoint& gpj);

  virtual void findIntersections(
      base::Grid& grid, size_t levelSum,
      std::unordered_map<size_t, std::shared_ptr<base::HashGridPoint>>& res);

  virtual void initializeCandidates(base::Grid& grid, std::vector<size_t>& negativeGridPoints);

  void computeIntersection(base::HashGridPoint& gpi, base::HashGridPoint& gpj,
                           base::HashGridPoint& gpintersection);

  static bool compareGridPointsByHash(const std::shared_ptr<base::HashGridPoint>& lhs,
                                      const std::shared_ptr<base::HashGridPoint>& rhs);

  std::unordered_map<size_t, std::shared_ptr<std::vector<std::shared_ptr<base::HashGridPoint>>>>
      intersections;
  std::unordered_map<size_t, std::shared_ptr<base::HashGridPoint>> currentIntersections;
  std::unordered_map<size_t, std::shared_ptr<base::HashGridPoint>> nextIntersections;
  std::unordered_map<size_t, std::shared_ptr<base::HashGridPoint>> candidates;

  size_t costs;
};

// -------------------------------------------------------------------------------------------
class OperationMakePositiveLoadFullGridCandidates
    : public OperationMakePositiveCandidateSetAlgorithm {
 public:
  OperationMakePositiveLoadFullGridCandidates();
  virtual ~OperationMakePositiveLoadFullGridCandidates();

  void nextCandidates(base::Grid& grid, base::DataVector& alpha, size_t levelSum,
                      std::vector<std::shared_ptr<base::HashGridPoint>>& candidates) override;

  size_t numCandidates() override;

 private:
  void initializeFullGrid(base::Grid& grid);

  std::unique_ptr<base::Grid> fullGrid;
};

// -------------------------------------------------------------------------------------------
class OperationMakePositiveHybridFindIntersectionCandidates
    : public OperationMakePositiveFindIntersectionCandidates {
 public:
  explicit OperationMakePositiveHybridFindIntersectionCandidates(size_t fullGridLevel);
  virtual ~OperationMakePositiveHybridFindIntersectionCandidates();

  void nextCandidates(base::Grid& grid, base::DataVector& alpha, size_t levelSum,
                      std::vector<std::shared_ptr<base::HashGridPoint>>& candidates) override;

 protected:
  void initializeCandidates(base::Grid& grid, std::vector<size_t>& negativeGridPoints) override;

  size_t fullGridLevel;
  std::unique_ptr<base::Grid> fullGrid;
};

} /* namespace datadriven */
} /* namespace sgpp */
