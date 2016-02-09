// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONMULTIPLEEVALPERIODIC_HPP
#define OPERATIONMULTIPLEEVALPERIODIC_HPP

#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {

/**
 * This class implements OperationMultipleEval for a grids with periodic linear basis ansatzfunctions
 *
 */
class OperationMultipleEvalPeriodic : public OperationMultipleEval {
 public:
  /**
   * Constructor
   *
   * @param grid grid
   * @param dataset the dataset that should be evaluated
   */
  OperationMultipleEvalPeriodic(Grid& grid,
                                DataMatrix& dataset) :
    OperationMultipleEval(grid, dataset) {
    this->storage = grid.getStorage();
  }

  /**
   * Destructor
   */
  ~OperationMultipleEvalPeriodic() override {}

  void mult(DataVector& alpha, DataVector& result) override;
  virtual void multTranspose(DataVector& source, DataVector& result);

 protected:
  /// Pointer to GridStorage object
  GridStorage* storage;
};

}  // namespace base
}  // namespace SGPP

#endif /* OPERATIONMULTIPLEEVALPERIODIC_HPP */
