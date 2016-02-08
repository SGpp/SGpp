// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONMULTIPLEEVALPREWAVELET_HPP
#define OPERATIONMULTIPLEEVALPREWAVELET_HPP

#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {

/**
 * This class implements OperationMultipleEval for a grids with prewavelet ansatzfunctions without boundaries
 */
class OperationMultipleEvalPrewavelet : public OperationMultipleEval {
 public:
  /**
   * Constructor of OperationMultipleEvalPrewavelet
   *
   * @param grid grid
   * @param dataset Dataset
   */
  OperationMultipleEvalPrewavelet(Grid& grid,
                                  DataMatrix& dataset) : OperationMultipleEval(grid, dataset) {
    this->storage = grid.getStorage();
  }

  /**
   * Destructor
   */
  virtual ~OperationMultipleEvalPrewavelet() override {}

  void mult(DataVector& alpha, DataVector& result) override;
  void multTranspose(DataVector& source, DataVector& result) override;

 protected:
  /// Pointer to the grid's GridStorage object
  GridStorage* storage;
};

}  // namespace base
}  // namespace SGPP

#endif /* OPERATIONMULTIPLEEVALPREWAVELET_HPP */
