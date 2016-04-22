// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONEVALPARTIALDERIVATIVEFUNDAMENTALSPLINE_HPP
#define OPERATIONEVALPARTIALDERIVATIVEFUNDAMENTALSPLINE_HPP

#include <sgpp/globaldef.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEvalPartialDerivative.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/operation/hash/common/basis/FundamentalSplineBasis.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

namespace sgpp {
namespace base {

/**
 * Operation for evaluating partial derivatives of B-spline
 * linear combinations on Noboundary grids.
 */
class OperationNaiveEvalPartialDerivativeFundamentalSpline :
  public OperationNaiveEvalPartialDerivative {
 public:
  /**
   * Constructor.
   *
   * @param storage   storage of the sparse grid
   * @param degree    B-spline degree
   */
  OperationNaiveEvalPartialDerivativeFundamentalSpline(GridStorage& storage,
      size_t degree) :
    storage(storage), base(degree) {
  }

  /**
   * Destructor.
   */
  ~OperationNaiveEvalPartialDerivativeFundamentalSpline() override {
  }

  /**
   * @param alpha     coefficient vector
   * @param point     evaluation point
   * @param derivDim  dimension in which the partial derivative should be taken
   * @return          value of the partial derivative of the linear combination
   */
  double evalPartialDerivative(const DataVector& alpha,
                                const DataVector& point,
                                size_t derivDim) override;

  /**
   * @param       alpha               coefficient matrix (each column is a coefficient vector)
   * @param       point               evaluation point
   * @param       derivDim            dimension in which the partial derivative should be taken
   *                                  (0, ..., d-1)
   * @param[out]  partialDerivative   values of the partial derivatives of the linear combination
   *                                  (the j-th entry corresponds to the j-th column of alpha)
   */
  void evalPartialDerivative(const DataMatrix& alpha,
                             const DataVector& point,
                             size_t derivDim,
                             DataVector& partialDerivative) override;

 protected:
  /// storage of the sparse grid
  GridStorage& storage;
  /// 1D B-spline basis
  SFundamentalSplineBase base;
};

}  // namespace base
}  // namespace sgpp

#endif /* OPERATIONEVALPARTIALDERIVATIVEFUNDAMENTALSPLINE_HPP */
