// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONQUADRATUREMODPOLY_HPP
#define OPERATIONQUADRATUREMODPOLY_HPP

#include <sgpp/base/operation/hash/OperationQuadrature.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/hash/common/basis/PolyModifiedBasis.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

/**
 * Quadrature on sparse grid, polynomial grid without boundaries
 */
class OperationQuadratureModPoly : public OperationQuadrature {
 public:
  /**
   * Constructor of OperationQuadraturePoly
   *
   * @param storage Pointer to the grid's GridStorage object
   * @param degree the polynom's max. degree
   */
  OperationQuadratureModPoly(GridStorage& storage,
                                  size_t degree) : storage(storage), base(degree) {}

  ~OperationQuadratureModPoly() override {}

  /**
   * Quadrature for Poly basis functions of max. degree 7
   *
   * @param alpha Coefficient vector for current grid
   */
  double doQuadrature(DataVector& alpha) override;

 protected:
  // Pointer to the grid's GridStorage object
  GridStorage& storage;
  /// Poly Boundary Basis object
  SPolyModifiedBase base;
};

}  // namespace base
}  // namespace sgpp

#endif /* OPERATIONQUADRATUREMODPOLY_HPP */
