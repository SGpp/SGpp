// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONQUADRATUREMODFUNDAMENTALSPLINE_HPP
#define OPERATIONQUADRATUREMODFUNDAMENTALSPLINE_HPP

#include <sgpp/base/operation/hash/OperationQuadrature.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/hash/common/basis/FundamentalSplineModifiedBasis.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

/**
 * Quadrature on sparse grid, ModFundamentalSplinenomial grid without boundaries
 */
class OperationQuadratureModFundamentalSpline : public OperationQuadrature {
 public:
  /**
   * Constructor of OperationQuadratureModFundamentalSpline
   *
   * @param storage Pointer to the grid's GridStorage object
   * @param degree the ModFundamentalSplines degree
   */
  OperationQuadratureModFundamentalSpline(GridStorage& storage, size_t degree)
      : storage(storage), base(degree) {}

  ~OperationQuadratureModFundamentalSpline() override {}

  /**
   * Quadrature for piecewise ModFundamentalSplinenomial basis functions of max. degree 3
   *
   * @param alpha Coefficient vector for current grid
   */
  double doQuadrature(DataVector& alpha) override;

 protected:
  // Pointer to the grid's GridStorage object
  GridStorage& storage;
  /// ModFundamentalSpline Basis object
  SFundamentalSplineModifiedBase base;
};

}  // namespace base
}  // namespace sgpp

#endif /* OPERATIONQUADRATUREMODFUNDAMENTALSPLINE_HPP */
