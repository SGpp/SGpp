// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/hash/OperationQuadrature.hpp>
#include <sgpp/globaldef.hpp>

#include "common/basis/NakPBsplineBasis.hpp"

namespace sgpp {
namespace base {

/**
 * Quadrature on sparse grid, not a knotter B-spline grid
 */
class OperationQuadratureNakPBspline : public OperationQuadrature {
 public:
  /**
   * Constructor of OperationQuadratureNakPBspline
   *
   * @param storage Pointer to the grid's GridStorage object
   * @param degree the B-spline degree
   */
  OperationQuadratureNakPBspline(GridStorage& storage, size_t degree)
      : storage(storage), base(degree) {}

  ~OperationQuadratureNakPBspline() override {}

  /**
   * Quadrature for not a knotter B-spline basis functions
   *
   * @param alpha Coefficient vector for current grid
   */
  double doQuadrature(DataVector& alpha) override;

 protected:
  // Pointer to the grid's GridStorage object
  GridStorage& storage;
  /// Basis object
  SNakPBsplineBase base;
};

}  // namespace base
}  // namespace sgpp
