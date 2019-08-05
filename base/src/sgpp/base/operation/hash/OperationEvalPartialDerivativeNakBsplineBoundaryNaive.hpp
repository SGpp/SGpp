// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONEVALPARTIALDERIVATIVENAKSPLINEBOUNDARY_HPP
#define OPERATIONEVALPARTIALDERIVATIVENAKSPLINEBOUNDARY_HPP

#include <sgpp/globaldef.hpp>
#include <sgpp/base/operation/hash/OperationEvalPartialDerivative.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineBasisDeriv1.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

namespace sgpp {
namespace base {

/**
 * Operation for evaluating partial derivatives of spline
 * linear combinations on Boundary grids with not-a-knot boundary conditions.
 */
class OperationEvalPartialDerivativeNakBsplineBoundaryNaive :
  public OperationEvalPartialDerivative {
 public:
  /**
   * Constructor.
   *
   * @param storage   storage of the sparse grid
   * @param degree    B-spline degree
   */
  OperationEvalPartialDerivativeNakBsplineBoundaryNaive(
      GridStorage& storage, size_t degree) :
    storage(storage),
    base(degree),
    baseDeriv1(degree),
    pointInUnitCube(storage.getDimension()) {
  }

  /**
   * Destructor.
   */
  ~OperationEvalPartialDerivativeNakBsplineBoundaryNaive() override {
  }

  /**
   * @param       alpha               coefficient vector
   * @param       point               evaluation point
   * @param       derivDim            dimension in which the partial derivative should be taken
   *                                  (0, ..., d-1)
   * @param[out]  partialDerivative   value of the partial derivative of the linear combination
   * @return                          value of the linear combination
   */
  double evalPartialDerivative(const DataVector& alpha,
                                const DataVector& point,
                                size_t derivDim,
                                double& partialDerivative) override;

  /**
   * @param       alpha               coefficient matrix (each column is a coefficient vector)
   * @param       point               evaluation point
   * @param       derivDim            dimension in which the partial derivative should be taken
   *                                  (0, ..., d-1)
   * @param[out]  value               values of the linear combination
   * @param[out]  partialDerivative   values of the partial derivatives of the linear combination
   *                                  (the j-th entry corresponds to the j-th column of alpha)
   */
  void evalPartialDerivative(const DataMatrix& alpha,
                             const DataVector& point,
                             size_t derivDim,
                             DataVector& value,
                             DataVector& partialDerivative) override;

 protected:
  /// storage of the sparse grid
  GridStorage& storage;
  /// 1D spline basis
  SNakBsplineBase base;
  /// 1D spline basis derivative
  SNakBsplineBaseDeriv1 baseDeriv1;
  /// untransformed evaluation point (temporary vector)
  DataVector pointInUnitCube;
};

}  // namespace base
}  // namespace sgpp

#endif /* OPERATIONEVALPARTIALDERIVATIVENAKSPLINEBOUNDARY_HPP */
