// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/operation/hash/OperationQuadrature.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/hash/common/basis/PolyClenshawCurtisBoundaryBasis.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

/**
 * Quadrature on sparse grid, polynomial grid without boundaries
 */
class OperationQuadraturePolyClenshawCurtisBoundary : public OperationQuadrature {
 public:
  /**
   * Constructor of OperationQuadraturePolyClenshawCurtisBoundary
   *
   * @param storage Pointer to the grid's GridStorage object
   * @param degree the polynom's max. degree
   */
  OperationQuadraturePolyClenshawCurtisBoundary(GridStorage& storage, size_t degree)
      : storage(storage), base(degree) {}

  ~OperationQuadraturePolyClenshawCurtisBoundary() override {}

  /**
   * Quadrature for piecewise polynomial basis functions of max. degree 3
   *
   * @param alpha Coefficient vector for current grid
   */
  double doQuadrature(DataVector& alpha) override;

 protected:
  // Pointer to the grid's GridStorage object
  GridStorage& storage;
  /// Poly Basis object
  SPolyClenshawCurtisBoundaryBase base;
};

}  // namespace base
}  // namespace sgpp
