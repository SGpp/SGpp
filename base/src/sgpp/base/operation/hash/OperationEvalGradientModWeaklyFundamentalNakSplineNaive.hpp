// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONEVALGRADIENTMODWEAKLYFUNDAMENTALNAKSPLINE_HPP
#define OPERATIONEVALGRADIENTMODWEAKLYFUNDAMENTALNAKSPLINE_HPP

#include <sgpp/globaldef.hpp>
#include <sgpp/base/operation/hash/OperationEvalGradient.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/operation/hash/common/basis/WeaklyFundamentalNakSplineModifiedBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/WeaklyFundamentalNakSplineModifiedBasisDeriv1.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

namespace sgpp {
namespace base {

/**
 * Operation for evaluating modified weakly fundamental spline linear combinations on Noboundary grids
 * and their gradients.
 */
class OperationEvalGradientModWeaklyFundamentalNakSplineNaive : public
  OperationEvalGradient {
 public:
  /**
   * Constructor.
   *
   * @param storage   storage of the sparse grid
   * @param degree    B-spline degree
   */
  OperationEvalGradientModWeaklyFundamentalNakSplineNaive(GridStorage& storage, size_t degree) :
    storage(storage),
    base(degree),
    baseDeriv1(degree),
    pointInUnitCube(storage.getDimension()),
    innerDerivative(storage.getDimension()) {
  }

  /**
   * Destructor.
   */
  ~OperationEvalGradientModWeaklyFundamentalNakSplineNaive() override {
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
  /// 1D spline basis
  SWeaklyFundamentalNakSplineModifiedBase base;
  /// 1D spline basis derivative
  SWeaklyFundamentalNakSplineModifiedBaseDeriv1 baseDeriv1;
  /// untransformed evaluation point (temporary vector)
  DataVector pointInUnitCube;
  /// inner derivative (temporary vector)
  DataVector innerDerivative;
};

}  // namespace base
}  // namespace sgpp

#endif /* OPERATIONEVALGRADIENTMODWEAKLYFUNDAMENTALNAKSPLINE_HPP */
