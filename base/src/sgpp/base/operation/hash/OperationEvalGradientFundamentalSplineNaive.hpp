// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONEVALGRADIENTFUNDAMENTALSPLINE_HPP
#define OPERATIONEVALGRADIENTFUNDAMENTALSPLINE_HPP

#include <sgpp/globaldef.hpp>
#include <sgpp/base/operation/hash/OperationEvalGradient.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/operation/hash/common/basis/FundamentalSplineBasis.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

namespace sgpp {
namespace base {

/**
 * Operation for evaluating fundamental sspline linear combinations on Noboundary grids and their
 * gradients.
 */
class OperationEvalGradientFundamentalSplineNaive : public
  OperationEvalGradient {
 public:
  /**
   * Constructor.
   *
   * @param storage   storage of the sparse grid
   * @param degree    B-spline degree
   */
  OperationEvalGradientFundamentalSplineNaive(GridStorage& storage, size_t degree) :
    storage(storage),
    base(degree),
    pointInUnitCube(storage.getDimension()),
    innerDerivative(storage.getDimension()) {
  }

  /**
   * Destructor.
   */
  ~OperationEvalGradientFundamentalSplineNaive() override {
  }

  /**
   * @param       alpha     coefficient vector
   * @param       point     evaluation point
   * @param[out]  gradient  gradient of linear combination
   * @return                value of the linear combination
   */
  double evalGradient(const DataVector& alpha,
                       const DataVector& point,
                       DataVector& gradient) override;

  /**
   * @param       alpha     coefficient matrix (each column is a coefficient vector)
   * @param       point     evaluation point
   * @param[out]  value     values of the linear combination
   * @param[out]  gradient  Jacobian of the linear combination (each row is a gradient vector)
   */
  void evalGradient(const DataMatrix& alpha,
                    const DataVector& point,
                    DataVector& value,
                    DataMatrix& gradient) override;

 protected:
  /// storage of the sparse grid
  GridStorage& storage;
  /// 1D B-spline basis
  SFundamentalSplineBase base;
  /// untransformed evaluation point (temporary vector)
  DataVector pointInUnitCube;
  /// inner derivative (temporary vector)
  DataVector innerDerivative;
};

}  // namespace base
}  // namespace sgpp

#endif /* OPERATIONEVALGRADIENTFUNDAMENTALSPLINE_HPP */
