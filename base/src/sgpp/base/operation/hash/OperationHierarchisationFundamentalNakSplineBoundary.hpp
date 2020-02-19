// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONHIERARCHISATIONFUNDAMENTALNAKSPLINEBOUNDARY_HPP
#define OPERATIONHIERARCHISATIONFUNDAMENTALNAKSPLINEBOUNDARY_HPP

#include <sgpp/base/operation/hash/OperationHierarchisation.hpp>
#include <sgpp/base/grid/type/FundamentalNakSplineBoundaryGrid.hpp>

#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace base {

/**
 * Hierarchisation on sparse grid, fundamental spline basis
 */
class OperationHierarchisationFundamentalNakSplineBoundary :
  public OperationHierarchisation {
 public:
  /**
   * Constructor of OperationHierarchisationFundamentalNakSplineBoundary
   *
   * @param grid Pointer to the grid
   */
  explicit OperationHierarchisationFundamentalNakSplineBoundary(
    FundamentalNakSplineBoundaryGrid* grid);

  /**
   * Destructor.
   */
  ~OperationHierarchisationFundamentalNakSplineBoundary() override;

  void doHierarchisation(DataVector& node_values) override;
  void doDehierarchisation(DataVector& alpha) override;

  void doHierarchisation(DataMatrix& node_values);
  void doDehierarchisation(DataMatrix& alpha);

 protected:
  /// grid
  FundamentalNakSplineBoundaryGrid* grid;
};

}  // namespace base
}  // namespace sgpp

#endif /* OPERATIONHIERARCHISATIONFUNDAMENTALNAKSPLINEBOUNDARY_HPP */
