// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OperationMultipleEvalMatrix_HPP
#define OperationMultipleEvalMatrix_HPP

#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>

namespace sgpp {
namespace datadriven {

/**
 * wrapper class for matrix multiplication for use in solver
 */
class OperationMultipleEvalMatrix : public base::OperationMultipleEval {
 public:
  /**
   * Constructor
   * @param grid dummy grid for inheritance reasons
   * @param dataset the dataset that should be evaluated
   */
  OperationMultipleEvalMatrix(base::Grid &grid, base::DataMatrix &dataset)
      : OperationMultipleEval(grid, dataset) {}

  /**
   * Destructor
   */
  ~OperationMultipleEvalMatrix() override {}

  void mult(base::DataVector &alpha, base::DataVector &result) override;
  void multTranspose(base::DataVector &source, base::DataVector &result) override;

  double getDuration() override;
};
}  // namespace datadriven
}  // namespace sgpp

#endif /* OperationMultipleEvalMatrix_HPP */
