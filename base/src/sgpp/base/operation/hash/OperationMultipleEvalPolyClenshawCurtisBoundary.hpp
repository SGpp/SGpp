// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp/base/operation/hash/common/basis/PolyClenshawCurtisBoundaryBasis.hpp>

namespace sgpp {
namespace base {

/**
 * This class implements OperationMultipleEval for a grids with poly basis ansatzfunctions
 */
class OperationMultipleEvalPolyClenshawCurtisBoundary : public OperationMultipleEval {
 public:
  /**
   * Constructor
   *
   * @param grid grid
   * @param degree the polynom's max. degree
   * @param dataset Dataset
   */
  OperationMultipleEvalPolyClenshawCurtisBoundary(Grid& grid, size_t degree, DataMatrix& dataset)
      : OperationMultipleEval(grid, dataset), storage(grid.getStorage()), base(degree) {}

  /**
   * Destructor
   */
  ~OperationMultipleEvalPolyClenshawCurtisBoundary() override {}

  void mult(DataVector& alpha, DataVector& result) override;
  void multTranspose(DataVector& source, DataVector& result) override;

  double getDuration() override;

 protected:
  /// Pointer to GridStorage object
  GridStorage& storage;
  /// Poly Basis object
  SPolyClenshawCurtisBoundaryBase base;
};

}  // namespace base
}  // namespace sgpp
