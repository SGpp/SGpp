// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONMULTIPLEEVALMODLINEAR_HPP
#define OPERATIONMULTIPLEEVALMODLINEAR_HPP

#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {

/**
 * This class implements OperationMultipleEval for a grids with mod linear basis ansatzfunctions
 *
 */
class OperationMultipleEvalModLinear : public OperationMultipleEval {
 public:
  /**
   * Constructor
   *
   * @param grid grid
   * @param dataset the dataset that should be evaluated
   */
  OperationMultipleEvalModLinear(Grid& grid,
                                 DataMatrix& dataset) :
                                   OperationMultipleEval(grid, dataset) {
    this->storage = grid.getStorage();
  }

  /**
   * Destructor
   */
  ~OperationMultipleEvalModLinear() override {}

  void mult(DataVector& alpha, DataVector& result) override;
  void multTranspose(DataVector& source, DataVector& result) override;

 protected:
  /// Pointer to GridStorage object
  GridStorage* storage;
};

}  // namespace base
}  // namespace SGPP

#endif /* OPERATIONMULTIPLEEVALMODLINEAR_HPP */
