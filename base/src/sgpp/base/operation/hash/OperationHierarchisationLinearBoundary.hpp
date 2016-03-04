// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONHIERARCHISATIONLINEARBOUNDARY_HPP
#define OPERATIONHIERARCHISATIONLINEARBOUNDARY_HPP

#include <sgpp/base/operation/hash/OperationHierarchisation.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace base {

/**
 * Hierarchisation on sparse grid, linear case with boundaries
 *
 */
class OperationHierarchisationLinearBoundary : public OperationHierarchisation {
 public:
  /**
   * Constructor
   *
   * @param storage the grid's GridStorage object
   */
  explicit OperationHierarchisationLinearBoundary(GridStorage& storage) :
    storage(storage) {}

  /**
   * Destructor
   */
  ~OperationHierarchisationLinearBoundary() override {}

  void doHierarchisation(DataVector& node_values) override;
  void doDehierarchisation(DataVector& alpha) override;

 protected:
  /// Pointer to GridStorage object
  GridStorage& storage;
};

}  // namespace base
}  // namespace sgpp

#endif /* OPERATIONHIERARCHISATIONLINEARBOUNDARY_HPP */
