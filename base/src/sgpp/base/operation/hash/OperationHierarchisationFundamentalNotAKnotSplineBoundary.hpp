// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONHIERARCHISATIONFUNDAMENTALNOTAKNOTSPLINEBOUNDARY_HPP
#define OPERATIONHIERARCHISATIONFUNDAMENTALNOTAKNOTSPLINEBOUNDARY_HPP

#include <sgpp/base/operation/hash/OperationHierarchisation.hpp>
#include <sgpp/base/grid/type/FundamentalNotAKnotSplineBoundaryGrid.hpp>

#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace base {

/**
 * Hierarchisation on sparse grid, fundamental spline basis
 */
class OperationHierarchisationFundamentalNotAKnotSplineBoundary :
  public OperationHierarchisation {
 public:
  /**
   * Constructor of OperationHierarchisationFundamentalNotAKnotSplineBoundary
   *
   * @param grid Pointer to the grid
   */
  explicit OperationHierarchisationFundamentalNotAKnotSplineBoundary(
    FundamentalNotAKnotSplineBoundaryGrid* grid);

  /**
   * Destructor.
   */
  ~OperationHierarchisationFundamentalNotAKnotSplineBoundary() override;

  void doHierarchisation(DataVector& node_values) override;
  void doDehierarchisation(DataVector& alpha) override;

  void doHierarchisation(DataMatrix& node_values);
  void doDehierarchisation(DataMatrix& alpha);

 protected:
  /// grid
  FundamentalNotAKnotSplineBoundaryGrid* grid;
};

}  // namespace base
}  // namespace sgpp

#endif /* OPERATIONHIERARCHISATIONFUNDAMENTALNOTAKNOTSPLINEBOUNDARY_HPP */
