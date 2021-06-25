// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/globaldef.hpp>
#include "common/basis/NakBsplineModifiedBasis.hpp"

namespace sgpp {
namespace base {

class OperationMultipleEvalNakBsplineModifiedNaive : public OperationMultipleEval {
 public:
  OperationMultipleEvalNakBsplineModifiedNaive(Grid& grid, size_t degree, DataMatrix& dataset)
      : OperationMultipleEval(grid, dataset), storage(grid.getStorage()), base(degree) {}

  ~OperationMultipleEvalNakBsplineModifiedNaive() override {}

  void mult(DataVector& alpha, DataVector& result) override;
  void multTranspose(DataVector& source, DataVector& result) override;

  double getDuration() override;

 protected:
  /// storage of the sparse grid
  GridStorage& storage;
  /// 1D B-spline basis
  SNakBsplineModifiedBase base;
  /// untransformed evaluation point (temporary vector)
  DataMatrix pointsInUnitCube;
};

}  // namespace base
}  // namespace sgpp
