// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONEVALGRADIENTMODFUNDAMENTALSPLINE_HPP
#define OPERATIONEVALGRADIENTMODFUNDAMENTALSPLINE_HPP

#include <sgpp/globaldef.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEvalGradient.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/operation/hash/common/basis/FundamentalSplineModifiedBasis.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

namespace SGPP {
namespace base {

/**
 * Operation for evaluating modified B-spline linear combinations on Noboundary grids and
 * their gradients.
 */
class OperationNaiveEvalGradientModFundamentalSpline : public
  OperationNaiveEvalGradient {
 public:
  /**
   * Constructor.
   *
   * @param storage   storage of the sparse grid
   * @param degree    B-spline degree
   */
  OperationNaiveEvalGradientModFundamentalSpline(GridStorage* storage,
      size_t degree) :
    storage(storage), base(degree) {
  }

  /**
   * Destructor.
   */
  ~OperationNaiveEvalGradientModFundamentalSpline() override {
  }

  /**
   * @param       alpha       coefficient vector
   * @param       point       evaluation point
   * @param[out]  gradient    gradient of linear combination
   * @return                  value of linear combination
   */
  float_t evalGradient(const DataVector& alpha,
                       const DataVector& point,
                       DataVector& gradient) override;

 protected:
  /// storage of the sparse grid
  GridStorage* storage;
  /// 1D B-spline basis
  SFundamentalSplineModifiedBase base;
};

}  // namespace base
}  // namespace SGPP

#endif /* OPERATIONEVALGRADIENTMODFUNDAMENTALSPLINE_HPP */
