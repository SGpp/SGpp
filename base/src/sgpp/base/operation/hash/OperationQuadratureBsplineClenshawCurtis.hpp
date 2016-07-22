// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONQUADRATUREBSPLINECLENSHAWCURTIS_HPP
#define OPERATIONQUADRATUREBSPLINECLENSHAWCURTIS_HPP

#include <sgpp/base/operation/hash/OperationQuadrature.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/hash/common/basis/BsplineClenshawCurtisBasis.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

/**
 * Quadrature on sparse grid, polynomial grid without boundaries
 */
class OperationQuadratureBsplineClenshawCurtis : public OperationQuadrature {
 public:
  /**
   * Constructor of OperationQuadratureBspline
   *
   * @param storage Pointer to the grid's GridStorage object
   * @param degree the polynom's max. degree
   */
  OperationQuadratureBsplineClenshawCurtis(GridStorage& storage,
                                  size_t degree) : storage(storage), base(degree) {}

  ~OperationQuadratureBsplineClenshawCurtis() override {}

  /**
   * Quadrature for Bspline basis functions of max. degree 7
   *
   * @param alpha Coefficient vector for current grid
   */
  double doQuadrature(DataVector& alpha) override;

 protected:
  // Pointer to the grid's GridStorage object
  GridStorage& storage;
  /// Bspline ClenshawCurtis Basis object
  SBsplineClenshawCurtisBase base;
};

}  // namespace base
}  // namespace sgpp

#endif /* OPERATIONQUADRATUREBSPLINECLENSHAWCURTIS_HPP */
