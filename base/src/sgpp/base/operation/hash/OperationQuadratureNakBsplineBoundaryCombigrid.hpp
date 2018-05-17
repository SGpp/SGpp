// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/hash/OperationQuadrature.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineBoundaryCombigridBasis.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

/**
 * Quadrature on sparse grid, not a knot B-spline grid created by transformation of a not a knot
 * B-spline combigrid
 */
class OperationQuadratureNakBsplineBoundaryCombigrid : public OperationQuadrature {
 public:
  /**
   * Constructor of OperationQuadratureNakBsplineBoundaryCombigrid
   *
   * @param storage Pointer to the grid's GridStorage object
   * @param degree the B-spline degree
   */
  OperationQuadratureNakBsplineBoundaryCombigrid(GridStorage& storage, size_t degree)
      : storage(storage), base(degree) {}

  ~OperationQuadratureNakBsplineBoundaryCombigrid() override {}

  /**
   * Quadrature for not a knot B-spline basis functions
   *
   * @param alpha Coefficient vector for current grid
   */
  double doQuadrature(DataVector& alpha) override;

 protected:
  // Pointer to the grid's GridStorage object
  GridStorage& storage;
  /// NakBsplineBoundaryCombigrid Basis object
  SNakBsplineBoundaryCombigridBase base;
};

}  // namespace base
}  // namespace sgpp
