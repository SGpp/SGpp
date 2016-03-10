// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONMULTIPLEEVALLINEARSTRETCHEDBOUNDARY_HPP
#define OPERATIONMULTIPLEEVALLINEARSTRETCHEDBOUNDARY_HPP

#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace base {

/**
 * This class implements OperationMultipleEval for a grids with linearstretched basis ansatzfunctions
 * with boundaries
 *
 */
class OperationMultipleEvalLinearStretchedBoundary : public
  OperationMultipleEval {
 public:
  /**
   * Constructor
   *
   * @param grid grid
   * @param dataset the dataset the should be evaluated
   */
  OperationMultipleEvalLinearStretchedBoundary(Grid& grid, DataMatrix& dataset) :
    OperationMultipleEval(grid, dataset), storage(grid.getStorage()) {
  }

  /**
   * Destructor
   */
  ~OperationMultipleEvalLinearStretchedBoundary() override {}

  void mult(DataVector& alpha, DataVector& result) override;
  void multTranspose(DataVector& source, DataVector& result) override;

 protected:
  /// Pointer to GridStorage object
  GridStorage& storage;
};

}  // namespace base
}  // namespace sgpp

#endif /* OPERATIONMULTIPLEEVALLINEARSTRETCHEDBOUNDARY_HPP */
