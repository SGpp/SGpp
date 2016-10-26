// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearClenshawCurtisBasis.hpp>

namespace sgpp {
namespace base {

/**
 * This class implements OperationMultipleEval for a grids with poly basis ansatzfunctions
 */
class OperationMultipleEvalLinearClenshawCurtis : public OperationMultipleEval {
 public:
  /**
   * Constructor
   *
   * @param grid grid
   * @param degree the polynom's max. degree
   * @param dataset Dataset
   */
  OperationMultipleEvalLinearClenshawCurtis(Grid& grid, DataMatrix& dataset)
      : OperationMultipleEval(grid, dataset), storage(grid.getStorage()) {}

  /**
   * Destructor
   */
  ~OperationMultipleEvalLinearClenshawCurtis() override {}

  void mult(DataVector& alpha, DataVector& result) override;
  void multTranspose(DataVector& source, DataVector& result) override;

  double getDuration() override;

 protected:
  /// Pointer to GridStorage object
  GridStorage& storage;
  /// Poly Basis object
  SLinearClenshawCurtisBase base;
};

}  // namespace base
}  // namespace sgpp
