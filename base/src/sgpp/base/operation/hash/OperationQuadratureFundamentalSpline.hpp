// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONQUADRATUREFUNDAMENTALSPLINE_HPP
#define OPERATIONQUADRATUREFUNDAMENTALSPLINE_HPP

#include <sgpp/base/operation/hash/OperationQuadrature.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/hash/common/basis/FundamentalSplineBasis.hpp>

#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace base {

/**
 * Quadrature on sparse grid, FundamentalSplinenomial grid without boundaries
 */
class OperationQuadratureFundamentalSpline : public OperationQuadrature {
 public:
  /**
   * Constructor of OperationQuadratureFundamentalSpline
   *
   * @param storage Pointer to the grid's GridStorage object
   * @param degree the FundamentalSplines degree
   */
  OperationQuadratureFundamentalSpline(GridStorage& storage, size_t degree) : storage(storage),
    base(degree) {}

  ~OperationQuadratureFundamentalSpline() override {}

  /**
   * Quadrature for piecewise FundamentalSplinenomial basis functions of max. degree 3
   *
   * @param alpha Coefficient vector for current grid
   */
  double doQuadrature(DataVector& alpha) override;

 protected:
  // Pointer to the grid's GridStorage object
  GridStorage& storage;
  /// FundamentalSpline Basis object
  SFundamentalSplineBase base;
};

}  // namespace base
}  // namespace sgpp

#endif /* OPERATIONQUADRATUREFUNDAMENTALSPLINE_HPP */
