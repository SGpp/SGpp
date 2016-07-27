// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONMULTIPLEEVALMODBSPLINENAIVE_HPP
#define OPERATIONMULTIPLEEVALMODBSPLINENAIVE_HPP

#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/operation/hash/common/basis/BsplineModifiedBasis.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

class OperationMultipleEvalModBsplineNaive: public OperationMultipleEval {
 public:
  OperationMultipleEvalModBsplineNaive(Grid& grid, size_t degree, DataMatrix& dataset) :
    OperationMultipleEval(grid, dataset), storage(grid.getStorage()), base(degree) {
  }

  ~OperationMultipleEvalModBsplineNaive() override {
  }

  void mult(DataVector& alpha, DataVector& result) override;
  void multTranspose(DataVector& source, DataVector& result) override;

 protected:
  /// storage of the sparse grid
  GridStorage& storage;
  /// 1D B-spline basis
  SBsplineModifiedBase base;
  /// untransformed evaluation point (temporary vector)
  DataMatrix pointsInUnitCube;
};

}  // namespace base
}  // namespace sgpp

#endif /* OPERATIONMULTIPLEEVALMODBSPLINENAIVE_HPP */
