// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONHIERARCHISATIONMODNOTAKNOTBSPLINE_HPP
#define OPERATIONHIERARCHISATIONMODNOTAKNOTBSPLINE_HPP

#include <sgpp/base/operation/hash/OperationHierarchisation.hpp>
#include <sgpp/base/grid/type/ModNotAKnotBsplineGrid.hpp>

#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace base {

/**
 * Hierarchisation on sparse grid, fundamental spline basis
 */
class OperationHierarchisationModNotAKnotBspline :
  public OperationHierarchisation {
 public:
  /**
   * Constructor of OperationHierarchisationModNotAKnotBspline
   *
   * @param grid Pointer to the grid
   */
  explicit OperationHierarchisationModNotAKnotBspline(
    ModNotAKnotBsplineGrid* grid);

  /**
   * Destructor.
   */
  ~OperationHierarchisationModNotAKnotBspline() override;

  void doHierarchisation(DataVector& node_values) override;
  void doDehierarchisation(DataVector& alpha) override;

  void doHierarchisation(DataMatrix& node_values);
  void doDehierarchisation(DataMatrix& alpha);

 protected:
  /// grid
  ModNotAKnotBsplineGrid* grid;
};

}  // namespace base
}  // namespace sgpp

#endif /* OPERATIONHIERARCHISATIONMODNOTAKNOTBSPLINE_HPP */
