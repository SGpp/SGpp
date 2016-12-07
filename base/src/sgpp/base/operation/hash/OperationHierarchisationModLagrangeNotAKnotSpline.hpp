// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONHIERARCHISATIONMODLAGRANGENOTAKNOTSPLINE_HPP
#define OPERATIONHIERARCHISATIONMODLAGRANGENOTAKNOTSPLINE_HPP

#include <sgpp/base/operation/hash/OperationHierarchisation.hpp>
#include <sgpp/base/grid/type/ModLagrangeNotAKnotSplineGrid.hpp>

#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace base {

/**
 * Hierarchisation on sparse grid, fundamental spline basis
 */
class OperationHierarchisationModLagrangeNotAKnotSpline :
  public OperationHierarchisation {
 public:
  /**
   * Constructor of OperationHierarchisationModLagrangeNotAKnotSpline
   *
   * @param grid Pointer to the grid
   */
  explicit OperationHierarchisationModLagrangeNotAKnotSpline(
      ModLagrangeNotAKnotSplineGrid* grid);

  /**
   * Destructor.
   */
  ~OperationHierarchisationModLagrangeNotAKnotSpline() override;

  void doHierarchisation(DataVector& node_values) override;
  void doDehierarchisation(DataVector& alpha) override;

  void doHierarchisation(DataMatrix& node_values);
  void doDehierarchisation(DataMatrix& alpha);

 protected:
  /// grid
  ModLagrangeNotAKnotSplineGrid* grid;
};

}  // namespace base
}  // namespace sgpp

#endif /* OPERATIONHIERARCHISATIONMODLAGRANGENOTAKNOTSPLINE_HPP */
