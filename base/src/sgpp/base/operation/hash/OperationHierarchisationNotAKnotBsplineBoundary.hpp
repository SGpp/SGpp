// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONHIERARCHISATIONNOTAKNOTBSPLINEBOUNDARY_HPP
#define OPERATIONHIERARCHISATIONNOTAKNOTBSPLINEBOUNDARY_HPP

#include <sgpp/base/operation/hash/OperationHierarchisation.hpp>
#include <sgpp/base/grid/type/NotAKnotBsplineBoundaryGrid.hpp>

#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace base {

/**
 * Hierarchisation on sparse grid, fundamental spline basis
 */
class OperationHierarchisationNotAKnotBsplineBoundary :
  public OperationHierarchisation {
 public:
  /**
   * Constructor of OperationHierarchisationNotAKnotBsplineBoundary
   *
   * @param grid Pointer to the grid
   */
  explicit OperationHierarchisationNotAKnotBsplineBoundary(
    NotAKnotBsplineBoundaryGrid* grid);

  /**
   * Destructor.
   */
  ~OperationHierarchisationNotAKnotBsplineBoundary() override;

  void doHierarchisation(DataVector& node_values) override;
  void doDehierarchisation(DataVector& alpha) override;

  void doHierarchisation(DataMatrix& node_values);
  void doDehierarchisation(DataMatrix& alpha);

 protected:
  /// grid
  NotAKnotBsplineBoundaryGrid* grid;
};

}  // namespace base
}  // namespace sgpp

#endif /* OPERATIONHIERARCHISATIONNOTAKNOTBSPLINEBOUNDARY_HPP */
