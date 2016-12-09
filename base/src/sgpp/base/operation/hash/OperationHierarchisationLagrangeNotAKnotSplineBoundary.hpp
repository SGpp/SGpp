// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONHIERARCHISATIONLAGRANGENOTAKNOTSPLINEBOUNDARY_HPP
#define OPERATIONHIERARCHISATIONLAGRANGENOTAKNOTSPLINEBOUNDARY_HPP

#include <sgpp/base/operation/hash/OperationHierarchisation.hpp>
#include <sgpp/base/grid/type/LagrangeNotAKnotSplineBoundaryGrid.hpp>

#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace base {

/**
 * Hierarchisation on sparse grid, Lagrange spline basis with not-a-knot boundary conditions
 */
class OperationHierarchisationLagrangeNotAKnotSplineBoundary :
  public OperationHierarchisation {
 public:
  /**
   * Constructor of OperationHierarchisationLagrangeNotAKnotSplineBoundary
   *
   * @param grid Pointer to the grid
   */
  explicit OperationHierarchisationLagrangeNotAKnotSplineBoundary(
    LagrangeNotAKnotSplineBoundaryGrid* grid);

  /**
   * Destructor.
   */
  ~OperationHierarchisationLagrangeNotAKnotSplineBoundary() override;

  void doHierarchisation(DataVector& node_values) override;
  void doDehierarchisation(DataVector& alpha) override;

  void doHierarchisation(DataMatrix& node_values);
  void doDehierarchisation(DataMatrix& alpha);

 protected:
  /// grid
  LagrangeNotAKnotSplineBoundaryGrid* grid;
};

}  // namespace base
}  // namespace sgpp

#endif /* OPERATIONHIERARCHISATIONLAGRANGENOTAKNOTSPLINEBOUNDARY_HPP */
