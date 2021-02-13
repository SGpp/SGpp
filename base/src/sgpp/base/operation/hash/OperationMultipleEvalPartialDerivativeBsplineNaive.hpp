/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * OperationMultipleEvalPartialDerivativeBsplineNaive.hpp
 *
 * Author: Paul Sarbu
 */

#pragma once

#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/operation/hash/common/basis/BsplineBasis.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

class OperationMultipleEvalPartialDerivativeBsplineNaive : public OperationMultipleEval {
 public:
  OperationMultipleEvalPartialDerivativeBsplineNaive(Grid& grid, size_t degree, DataMatrix& dataset,
                                                     size_t derivD)
      : OperationMultipleEval(grid, dataset),
        storage(grid.getStorage()),
        base(degree),
        derivDim(derivD) {}

  ~OperationMultipleEvalPartialDerivativeBsplineNaive() override {}

  void mult(DataVector& alpha, DataVector& result) override;
  void multTranspose(DataVector& source, DataVector& result) override;

  double getDuration() override;

 protected:
  /// storage of the sparse grid
  GridStorage& storage;
  /// 1D B-spline basis
  SBsplineBase base;
  /// untransformed evaluation point (temporary vector)
  DataMatrix pointsInUnitCube;
  /// dimension of partial derivative
  size_t derivDim;
};

}  // namespace base
}  // namespace sgpp
