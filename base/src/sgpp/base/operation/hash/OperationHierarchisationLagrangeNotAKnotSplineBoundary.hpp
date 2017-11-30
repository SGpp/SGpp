// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONHIERARCHISATIONLAGRANGENOTAKNOTSPLINEBOUNDARY_HPP
#define OPERATIONHIERARCHISATIONLAGRANGENOTAKNOTSPLINEBOUNDARY_HPP

#include <sgpp/base/operation/hash/OperationHierarchisation.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace base {

/**
 * Hierarchisation on sparse grid, Lagrange not-a-knot spline case with boundaries
 *
 */
class OperationHierarchisationLagrangeNotAKnotSplineBoundary : public OperationHierarchisation {
 public:
  /**
   * Constructor
   *
   * @param storage the grid's GridStorage object
   */
  explicit OperationHierarchisationLagrangeNotAKnotSplineBoundary(GridStorage& storage) :
    storage(storage) {}

  /**
   * Destructor
   */
  ~OperationHierarchisationLagrangeNotAKnotSplineBoundary() override {}

  void doHierarchisation(DataVector& node_values) override;
  void doDehierarchisation(DataVector& alpha) override;

 protected:
  /// Pointer to GridStorage object
  GridStorage& storage;
};

}  // namespace base
}  // namespace sgpp

#endif /* OPERATIONHIERARCHISATIONLAGRANGENOTAKNOTSPLINEBOUNDARY_HPP */
