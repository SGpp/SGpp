// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONHIERARCHISATIONLAGRANGESPLINEBOUNDARY_HPP
#define OPERATIONHIERARCHISATIONLAGRANGESPLINEBOUNDARY_HPP

#include <sgpp/base/operation/hash/OperationHierarchisation.hpp>
#include <sgpp/base/grid/type/LagrangeSplineBoundaryGrid.hpp>

#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace base {

/**
 * Hierarchisation on sparse grid, fundamental spline basis
 */
class OperationHierarchisationLagrangeSplineBoundary :
  public OperationHierarchisation {
 public:
  /**
   * Constructor of OperationHierarchisationLagrangeSplineBoundary
   *
   * @param grid Pointer to the grid
   */
  explicit OperationHierarchisationLagrangeSplineBoundary(
    LagrangeSplineBoundaryGrid* grid);

  /**
   * Destructor.
   */
  ~OperationHierarchisationLagrangeSplineBoundary() override;

  void doHierarchisation(DataVector& node_values) override;
  void doDehierarchisation(DataVector& alpha) override;

  void doHierarchisation(DataMatrix& node_values);
  void doDehierarchisation(DataMatrix& alpha);

 protected:
  /// grid
  LagrangeSplineBoundaryGrid* grid;
};

}  // namespace base
}  // namespace sgpp

#endif /* OPERATIONHIERARCHISATIONLAGRANGESPLINEBOUNDARY_HPP */
