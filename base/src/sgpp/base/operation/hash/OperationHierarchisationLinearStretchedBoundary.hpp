// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONHIERARCHISATIONLINEARSTRETCHEDBOUNDARY_HPP
#define OPERATIONHIERARCHISATIONLINEARSTRETCHEDBOUNDARY_HPP

#include <sgpp/base/operation/hash/OperationHierarchisation.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {

/**
 * Hierarchisation on sparse grid, linear stretched case with boundaries
 *
 */
class OperationHierarchisationLinearStretchedBoundary : public
  OperationHierarchisation {
 public:
  /**
   * Constructor
   *
   * @param storage the grid's GridStorage object
   */
  explicit OperationHierarchisationLinearStretchedBoundary(
    GridStorage& storage) : storage(storage) {}

  /**
   * Destructor
   */
  ~OperationHierarchisationLinearStretchedBoundary() override {}

  void doHierarchisation(DataVector& node_values) override;
  void doDehierarchisation(DataVector& alpha) override;

 protected:
  /// Pointer to GridStorage object
  GridStorage& storage;
};

}  // namespace base
}  // namespace SGPP

#endif /* OPERATIONHIERARCHISATIONLINEARSTRETCHEDBOUNDARY_HPP */
