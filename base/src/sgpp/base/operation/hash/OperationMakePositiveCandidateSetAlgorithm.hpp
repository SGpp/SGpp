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
namespace base {

class OperationMakePositiveCandidateSetAlgorithm {
 public:
  OperationMakePositiveCandidateSetAlgorithm();
  virtual ~OperationMakePositiveCandidateSetAlgorithm();

  virtual void nextCandidates(base::Grid& grid, base::DataVector& alpha, size_t levelSum,
                              std::vector<std::shared_ptr<HashGridPoint>>& candidates) = 0;

 protected:
  void findNodesWithNegativeCoefficients(base::DataVector& alpha,
                                         std::vector<size_t>& negativeGridPoints);

  bool hasOverlappingSupport(HashGridPoint& gpi, HashGridPoint& gpj, size_t dim);
  bool hasOverlappingSupport(HashGridPoint& gpi, HashGridPoint& gpj);
  bool isHierarchicalAncestor(HashGridPoint& gpi, HashGridPoint& gpj);
  bool isHierarchicalAncestor(HashGridPoint& gpi, HashGridPoint& gpj, size_t dim);
};
// -------------------------------------------------------------------------------------------
class OperationMakePositiveFindIntersectionCandidates
    : public OperationMakePositiveCandidateSetAlgorithm {
 public:
  explicit OperationMakePositiveFindIntersectionCandidates(base::Grid& grid);
  virtual ~OperationMakePositiveFindIntersectionCandidates();

  void nextCandidates(base::Grid& grid, base::DataVector& alpha, size_t levelSum,
                      std::vector<std::shared_ptr<HashGridPoint>>& candidates) override;

 private:
  void findIntersections(base::Grid& grid, size_t levelSum,
                         std::unordered_map<size_t, std::shared_ptr<HashGridPoint>>& res);

  void initializeCandidates(base::Grid& grid, std::vector<size_t>& negativeGridPoints);

  void computeIntersection(base::HashGridPoint& gpi, base::HashGridPoint& gpj,
                           base::HashGridPoint& gpintersection);

  size_t iteration;
  std::unordered_map<size_t, std::shared_ptr<std::vector<std::shared_ptr<HashGridPoint>>>>
      intersections;
  std::unordered_map<size_t, std::shared_ptr<HashGridPoint>> currentIntersections;
  std::unordered_map<size_t, std::shared_ptr<HashGridPoint>> nextIntersections;
  std::unordered_map<size_t, std::shared_ptr<HashGridPoint>> candidates;

  size_t costs;
  bool verbose;
};
// -------------------------------------------------------------------------------------------
class OperationMakePositiveLoadFullGridCandidates
    : public OperationMakePositiveCandidateSetAlgorithm {
 public:
  explicit OperationMakePositiveLoadFullGridCandidates(base::Grid& grid);
  virtual ~OperationMakePositiveLoadFullGridCandidates();

  void nextCandidates(base::Grid& grid, base::DataVector& alpha, size_t levelSum,
                      std::vector<std::shared_ptr<HashGridPoint>>& candidates) override;

 private:
  std::unique_ptr<base::Grid> fullGrid;
};

} /* namespace base */
} /* namespace sgpp */
