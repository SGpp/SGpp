// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONQUADRATUREBSPLINE_HPP
#define OPERATIONQUADRATUREBSPLINE_HPP

#include <sgpp/base/operation/hash/OperationQuadrature.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/hash/common/basis/BsplineBasis.hpp>

#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace base {

/**
 * Quadrature on sparse grid, Bsplinenomial grid without boundaries
 */
class OperationQuadratureBspline : public OperationQuadrature {
 public:
  /**
   * Constructor of OperationQuadratureBspline
   *
   * @param storage Pointer to the grid's GridStorage object
   * @param degree the Bsplines degree
   */
  OperationQuadratureBspline(GridStorage& storage, size_t degree) : storage(storage),
    base(degree) {}

  ~OperationQuadratureBspline() override {}

  /**
   * Quadrature for piecewise Bsplinenomial basis functions of max. degree 3
   *
   * @param alpha Coefficient vector for current grid
   */
  double doQuadrature(DataVector& alpha) override;

 protected:
  // Pointer to the grid's GridStorage object
  GridStorage& storage;
  /// Bspline Basis object
  SBsplineBase base;
};

}  // namespace base
}  // namespace sgpp

#endif /* OPERATIONQUADRATUREBSPLINE_HPP */
