// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/hash/OperationQuadrature.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineBoundaryBasis.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

/**
 * Quadrature on sparse grid, not a knot boundary B-spline grid
 */
class OperationQuadratureNakBsplineBoundary : public OperationQuadrature {
 public:
  /**
   * Constructor of OperationQuadratureNakBsplineBoundary
   *
   * @param storage Pointer to the grid's GridStorage object
   * @param degree the B-spline degree
   */
  OperationQuadratureNakBsplineBoundary(GridStorage& storage, size_t degree)
      : storage(storage), base(degree) {}

  ~OperationQuadratureNakBsplineBoundary() override {}

  /**
   * Quadrature for not a knot B-spline basis functions
   *
   * @param alpha Coefficient vector for current grid
   */
  double doQuadrature(DataVector& alpha) override;

 protected:
  // Pointer to the grid's GridStorage object
  GridStorage& storage;
  ///  Basis object
  SNakBsplineBoundaryBase base;
};

}  // namespace base
}  // namespace sgpp
