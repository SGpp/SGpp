// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONMULTIPLEEVALPOLYBOUNDARY_HPP
#define OPERATIONMULTIPLEEVALPOLYBOUNDARY_HPP

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/operation/hash/common/basis/PolyBoundaryBasis.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

/**
 * This class implements OperationMultipleEval for a grids with poly basis ansatzfunctions
 */
class OperationMultipleEvalPolyBoundary : public OperationMultipleEval {
 public:
  /**
   * Constructor
   *
   * @param grid the grid
   * @param degree the polynom's max. degree
   * @param dataset Dataset
   */
  OperationMultipleEvalPolyBoundary(Grid& grid, size_t degree, DataMatrix& dataset)
      : OperationMultipleEval(grid, dataset), storage(grid.getStorage()), base(degree) {}

  /**
   * Destructor
   */
  ~OperationMultipleEvalPolyBoundary() override {}

  void mult(DataVector& alpha, DataVector& result) override;
  void multTranspose(DataVector& source, DataVector& result) override;

  double getDuration() override;

 protected:
  /// Pointer to GridStorage object
  GridStorage& storage;
  /// Poly Basis object
  SPolyBoundaryBase base;
};

}  // namespace base
}  // namespace sgpp

#endif /* OPERATIONMULTIPLEEVALPOLYBOUNDARY_HPP */
