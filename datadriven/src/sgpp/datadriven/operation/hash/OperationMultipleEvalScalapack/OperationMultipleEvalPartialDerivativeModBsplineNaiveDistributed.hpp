/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * OperationMultipleEvalPartialDerivativeModBsplineNaive.hpp
 *
 */

#pragma once

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/operation/hash/common/basis/BsplineModifiedBasis.hpp>
#include <sgpp/datadriven/operation/hash/OperationMultipleEvalScalapack/OperationMultipleEvalDistributed.hpp>
#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {

class OperationMultipleEvalPartialDerivativeModBsplineNaiveDistributed
    : public OperationMultipleEvalDistributed {
 public:
  /**
   * Constructor
   * @param grid
   * @param degree
   * @param dataset
   * @param derivD
   */
  OperationMultipleEvalPartialDerivativeModBsplineNaiveDistributed(base::Grid& grid, size_t degree,
                                                                   DataMatrix& dataset,
                                                                   size_t derivD)
      : OperationMultipleEvalDistributed(grid, dataset),
        storage(grid.getStorage()),
        base(degree),
        derivDim(derivD) {}

  ~OperationMultipleEvalPartialDerivativeModBsplineNaiveDistributed() override {}

  /**
   * Distributed version of mult.
   *
   * @param alpha vector, to which @f$B@f$ is applied. Typically the coefficient vector
   * @param result the result vector of the matrix vector multiplication
   */
  void multDistributed(DataVector& alpha, DataVectorDistributed& result) override;
  /**
   * Distributed version of multTransposed.
   * Multiplication of @f$B@f$ with vector @f$\alpha@f$
   *
   * @param source vector, to which @f$B^T@f$ is applied. Typically the coefficient vector
   * @param result the result vector of the matrix vector multiplication
   */
  void multTransposeDistributed(DataVector& source, DataVectorDistributed& result) override;

  double getDuration() override;

 protected:
  /// storage of the sparse grid
  base::GridStorage& storage;
  /// 1D B-spline basis
  base::SBsplineModifiedBase base;
  /// untransformed evaluation point (temporary vector)
  DataMatrix pointsInUnitCube;
  /// dimension of partial derivative
  size_t derivDim;
};

}  // namespace datadriven
}  // namespace sgpp
