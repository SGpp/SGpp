// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/globaldef.hpp>

#include <vector>

namespace sgpp {
namespace base {

class OperationMakePositiveCandidateSetAlgorithm {
 public:
  OperationMakePositiveCandidateSetAlgorithm();
  virtual ~OperationMakePositiveCandidateSetAlgorithm();

  virtual void nextCandidates(base::Grid& grid, base::DataVector& alpha, size_t levelSum,
                              std::vector<HashGridPoint*>& candidates) = 0;
};
// -------------------------------------------------------------------------------------------
class OperationMakePositiveFindIntersectionCandidates
    : public OperationMakePositiveCandidateSetAlgorithm {
 public:
  OperationMakePositiveFindIntersectionCandidates();
  virtual ~OperationMakePositiveFindIntersectionCandidates();

  void nextCandidates(base::Grid& grid, base::DataVector& alpha, size_t levelSum,
                      std::vector<HashGridPoint*>& candidates) override;
};
// -------------------------------------------------------------------------------------------
class OperationMakePositiveLoadFullGridCandidates
    : public OperationMakePositiveCandidateSetAlgorithm {
 public:
  explicit OperationMakePositiveLoadFullGridCandidates(base::Grid& grid);
  virtual ~OperationMakePositiveLoadFullGridCandidates();

  void nextCandidates(base::Grid& grid, base::DataVector& alpha, size_t levelSum,
                      std::vector<HashGridPoint*>& candidates) override;

 private:
  std::unique_ptr<base::Grid> fullGrid;
};

} /* namespace base */
} /* namespace sgpp */
