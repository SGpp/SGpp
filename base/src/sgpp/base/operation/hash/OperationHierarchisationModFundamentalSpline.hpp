// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONHIERARCHISATIONMODFUNDAMENTALSPLINE_HPP
#define OPERATIONHIERARCHISATIONMODFUNDAMENTALSPLINE_HPP

#include <sgpp/base/operation/hash/OperationHierarchisation.hpp>
#include <sgpp/base/grid/type/ModFundamentalSplineGrid.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {

/**
 * Hierarchisation on sparse grid, modified fundamental spline basis
 */
class OperationHierarchisationModFundamentalSpline :
  public OperationHierarchisation {
 public:
  /**
   * Constructor of OperationHierarchisationModFundamentalSpline
   *
   * @param grid Pointer to the grid
   */
  OperationHierarchisationModFundamentalSpline(
    ModFundamentalSplineGrid* grid);

  /**
   * Destructor.
   */
  ~OperationHierarchisationModFundamentalSpline() override;

  void doHierarchisation(DataVector& node_values) override;
  void doDehierarchisation(DataVector& alpha) override;

  void doHierarchisation(DataMatrix& node_values);
  void doDehierarchisation(DataMatrix& alpha);

 protected:
  /// grid
  ModFundamentalSplineGrid* grid;
};

}  // namespace base
}  // namespace SGPP

#endif /* OPERATIONHIERARCHISATIONMODFUNDAMENTALSPLINE_HPP */
