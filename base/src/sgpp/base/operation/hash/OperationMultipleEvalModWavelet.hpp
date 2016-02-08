// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONMULTIPLEEVALMODWAVELET_HPP
#define OPERATIONMULTIPLEEVALMODWAVELET_HPP

#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {

/**
 * This class implements OperationMultipleEval for a grids with mod wavelet basis ansatzfunctions
 */
class OperationMultipleEvalModWavelet : public OperationMultipleEval {
 public:
  /**
   * Constructor
   *
   * @param grid grid
   * @param dataset data
   */
  OperationMultipleEvalModWavelet(Grid& grid,
                                  DataMatrix& dataset) : OperationMultipleEval(grid, dataset) {
    this->storage = grid.getStorage();
  }

  /**
   * Destructor
   */
  virtual ~OperationMultipleEvalModWavelet() override {}

  void mult(DataVector& alpha, DataVector& result) override;
  void multTranspose(DataVector& source, DataVector& result) override;

 protected:
  /// Pointer to GridStorage object
  GridStorage* storage;
};

}  // namespace base
}  // namespace SGPP

#endif /* OPERATIONMULTIPLEEVALMODWAVELET_HPP */
