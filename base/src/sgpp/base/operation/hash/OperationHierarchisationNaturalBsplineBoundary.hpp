// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONHIERARCHISATIONNATURALBSPLINEBOUNDARY_HPP
#define OPERATIONHIERARCHISATIONNATURALBSPLINEBOUNDARY_HPP

#include <sgpp/base/operation/hash/OperationHierarchisation.hpp>
#include <sgpp/base/grid/type/NaturalBsplineBoundaryGrid.hpp>

#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace base {

/**
 * Hierarchisation on sparse grid, fundamental spline basis
 */
class OperationHierarchisationNaturalBsplineBoundary :
  public OperationHierarchisation {
 public:
  /**
   * Constructor of OperationHierarchisationNaturalBsplineBoundary
   *
   * @param grid Pointer to the grid
   */
  explicit OperationHierarchisationNaturalBsplineBoundary(
    NaturalBsplineBoundaryGrid* grid);

  /**
   * Destructor.
   */
  ~OperationHierarchisationNaturalBsplineBoundary() override;

  void doHierarchisation(DataVector& node_values) override;
  void doDehierarchisation(DataVector& alpha) override;

  void doHierarchisation(DataMatrix& node_values);
  void doDehierarchisation(DataMatrix& alpha);

 protected:
  /// grid
  NaturalBsplineBoundaryGrid* grid;
};

}  // namespace base
}  // namespace sgpp

#endif /* OPERATIONHIERARCHISATIONNATURALBSPLINEBOUNDARY_HPP */
