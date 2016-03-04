// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONHIERARCHISATIONLINEARSTRETCHED_HPP
#define OPERATIONHIERARCHISATIONLINEARSTRETCHED_HPP

#include <sgpp/base/operation/hash/OperationHierarchisation.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace base {

/**
 * Hierarchisation on sparse grid, linear grid without boundaries
 */
class OperationHierarchisationLinearStretched : public
  OperationHierarchisation {
 public:
  /**
   * Constructor of OperationHierarchisationLinear
   *
   * @param storage Pointer to the grid's gridstorage obejct
   */
  explicit OperationHierarchisationLinearStretched(GridStorage& storage) :
    storage(storage) {}

  /**
   * Destructor
   */
  ~OperationHierarchisationLinearStretched() override {}

  void doHierarchisation(DataVector& node_values) override;
  void doDehierarchisation(DataVector& alpha) override;

 protected:
  /// reference to the grid's GridStorage object
  GridStorage& storage;
};

}  // namespace base
}  // namespace sgpp

#endif /* OPERATIONHIERARCHISATIONLINEARSTRETCHED_HPP */
