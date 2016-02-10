// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONQUADRATUREPOLYBOUNDARY_HPP
#define OPERATIONQUADRATUREPOLYBOUNDARY_HPP

#include <sgpp/base/operation/hash/OperationQuadrature.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/hash/common/basis/PolyBoundaryBasis.hpp>

#include <sgpp/globaldef.hpp>

namespace SGPP {
namespace base {

/**
 * Quadrature on sparse grid, polynomial grid without boundaries
 */
class OperationQuadraturePolyBoundary : public OperationQuadrature {
 public:
  /**
   * Constructor of OperationQuadraturePoly
   *
   * @param storage Pointer to the grid's GridStorage object
   * @param degree the polynom's max. degree
   */
  OperationQuadraturePolyBoundary(GridStorage* storage,
                                  size_t degree) : storage(storage), base(degree) {}

  ~OperationQuadraturePolyBoundary() override {}

  /**
   * Quadrature for piecewise polynomial basis functions of max. degree 3
   *
   * @param alpha Coefficient vector for current grid
   */
  float_t doQuadrature(DataVector& alpha) override;

 protected:
  // Pointer to the grid's GridStorage object
  GridStorage* storage;
  /// Poly Boundary Basis object
  SPolyBoundaryBase base;
};

}  // namespace base
}  // namespace SGPP

#endif /* OPERATIONQUADRATUREPOLYBOUNDARY_HPP */
