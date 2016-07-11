// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include "OperationMakePositiveCandidateSetAlgorithm.hpp"
#include <sgpp/base/operation/BaseOpFactory.hpp>

#include <vector>

namespace sgpp {
namespace base {

OperationMakePositiveCandidateSetAlgorithm::OperationMakePositiveCandidateSetAlgorithm() {}

OperationMakePositiveCandidateSetAlgorithm::~OperationMakePositiveCandidateSetAlgorithm() {}

// -------------------------------------------------------------------------------------------
OperationMakePositiveFindIntersectionCandidates::OperationMakePositiveFindIntersectionCandidates() {
}
OperationMakePositiveFindIntersectionCandidates::
    ~OperationMakePositiveFindIntersectionCandidates() {}

void OperationMakePositiveFindIntersectionCandidates::nextCandidates(
    base::Grid& grid, base::DataVector& alpha, size_t levelSum,
    std::vector<HashGridPoint*>& candidates) {}
// -------------------------------------------------------------------------------------------
OperationMakePositiveLoadFullGridCandidates::OperationMakePositiveLoadFullGridCandidates(
    base::Grid& grid) {
  base::HashGridStorage& gridStorage = grid.getStorage();

  fullGrid = Grid::createLinearGrid(gridStorage.getDimension());
  fullGrid->getGenerator().full(gridStorage.getMaxLevel());
}

OperationMakePositiveLoadFullGridCandidates::~OperationMakePositiveLoadFullGridCandidates() {}

void OperationMakePositiveLoadFullGridCandidates::nextCandidates(
    base::Grid& grid, base::DataVector& alpha, size_t levelSum,
    std::vector<HashGridPoint*>& candidates) {
  base::HashGridStorage& fullGridStorage = fullGrid->getStorage();
  base::HashGridPoint gp;
  for (size_t i = 0; i < fullGridStorage.getSize(); ++i) {
    gp = fullGridStorage.getPoint(i);
    if (gp.getLevelSum() == levelSum) {
      candidates.push_back(&gp);
    }
  }
}

} /* namespace base */
} /* namespace sgpp */
