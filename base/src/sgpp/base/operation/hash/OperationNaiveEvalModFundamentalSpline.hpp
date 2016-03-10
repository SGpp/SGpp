// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONNAIVEEVALMODFUNDAMENTALSPLINE_HPP
#define OPERATIONNAIVEEVALMODFUNDAMENTALSPLINE_HPP

#include <sgpp/globaldef.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEval.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/operation/hash/common/basis/FundamentalSplineModifiedBasis.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

namespace sgpp {
namespace base {

/**
 * Operation for evaluating modified B-spline linear combinations on Noboundary grids.
 */
class OperationNaiveEvalModFundamentalSpline : public OperationNaiveEval {
 public:
  /**
   * Constructor.
   *
   * @param storage   storage of the sparse grid
   * @param degree    B-spline degree
   */
  OperationNaiveEvalModFundamentalSpline(GridStorage& storage, size_t degree) :
    storage(storage), base(degree) {
  }

  /**
   * Destructor.
   */
  ~OperationNaiveEvalModFundamentalSpline() override {
  }

  /**
   * @param alpha     coefficient vector
   * @param point     evaluation point
   * @return          value of linear combination
   */
  double eval(const DataVector& alpha, const DataVector& point) override;

 protected:
  /// storage of the sparse grid
  GridStorage& storage;
  /// 1D B-spline basis
  SFundamentalSplineModifiedBase base;
};

}  // namespace base
}  // namespace sgpp

#endif /* OPERATIONNAIVEEVALMODFUNDAMENTALSPLINE_HPP */
