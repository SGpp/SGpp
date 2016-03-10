// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONMULTIPLEEVALLINEAR_HPP
#define OPERATIONMULTIPLEEVALLINEAR_HPP

#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace base {

/**
 * This class implements OperationB for a grids with linear basis ansatzfunctions without boundaries
 */
class OperationMultipleEvalLinear: public OperationMultipleEval {
 public:
  /**
   * Constructor of OperationBLinear
   *
   * @param grid grid
   * @param dataset the dataset that should be evaluated
   */
  OperationMultipleEvalLinear(Grid& grid, DataMatrix& dataset) :
    OperationMultipleEval(grid, dataset), storage(grid.getStorage()) {
  }

  /**
   * Destructor
   */
  ~OperationMultipleEvalLinear() override {
  }

  void mult(DataVector& alpha, DataVector& result) override;
  void multTranspose(DataVector& source, DataVector& result) override;

 protected:
  /// reference to the grid's GridStorage object
  GridStorage& storage;
};

}  // namespace base
}  // namespace sgpp

#endif /* OPERATIONMULTIPLEEVALLINEAR_HPP */
