/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * OperationMultipleEvalDistributed.hpp
 *
 * Created on: Mar 25, 2019
 *     Author: Jan Schopohl
 */

#pragma once

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/exception/not_implemented_exception.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/datadriven/scalapack/DataVectorDistributed.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {

/**
 * Distributed interface for OperationMultipleEval
 */
class OperationMultipleEvalDistributed : public sgpp::base::OperationMultipleEval {
 public:
  /**
   * Constructor
   *
   * @param grid the sparse grid used for this operation
   * @param dataset dataset that should be evaluated on the sparse grid
   */
  OperationMultipleEvalDistributed(sgpp::base::Grid& grid, sgpp::base::DataMatrix& dataset)
      : sgpp::base::OperationMultipleEval(grid, dataset) {}

  /**
   * Destructor
   */
  ~OperationMultipleEvalDistributed() override {}

  /**
   * Distributed version of mult.
   * Multiplication of @f$B^T@f$ with vector @f$\alpha@f$
   *
   * @param alpha vector, to which @f$B@f$ is applied. Typically the coefficient vector
   * @param result the result vector of the matrix vector multiplication
   */
  virtual void multDistributed(sgpp::base::DataVector& alpha, DataVectorDistributed& result) = 0;

  /**
   * Distributed version of multTransposed.
   * Multiplication of @f$B@f$ with vector @f$\alpha@f$
   *
   * @param source vector, to which @f$B^T@f$ is applied. Typically the coefficient vector
   * @param result the result vector of the matrix vector multiplication
   */
  virtual void multTransposeDistributed(sgpp::base::DataVector& source,
                                        DataVectorDistributed& result) = 0;

  /**
   * Evaluate multiple datapoints with the specified grid in parallel.
   *
   * @param alpha surplus vector of the grid
   * @param result distributed result of the evaluations
   */
  void evalParallel(DataVector& alpha, DataVectorDistributed& result) {
    this->multDistributed(alpha, result);
  }

  /**
   * Not implemented in distributed version.
   */
  virtual void mult(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result) {
    throw sgpp::base::not_implemented_exception("only distributed version of mult implemented");
  }

  /**
   * Not implemented in distributed version.
   */
  virtual void multTranspose(sgpp::base::DataVector& source, sgpp::base::DataVector& result) {
    throw sgpp::base::not_implemented_exception("only distributed version of mult implemented");
  }
};

}  // namespace datadriven
}  // namespace sgpp
