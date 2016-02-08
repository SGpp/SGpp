// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONHIERARCHISATIONLINEAR_HPP
#define OPERATIONHIERARCHISATIONLINEAR_HPP

#include <sgpp/base/operation/hash/OperationHierarchisation.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {

/**
 * Hierarchisation on sparse grid, linear grid without boundaries
 */
class OperationHierarchisationLinear : public OperationHierarchisation {
 public:
  /**
   * Constructor of OperationHierarchisationLinear
   *
   * @param storage Pointer to the grid's gridstorage obejct
   */
  explicit OperationHierarchisationLinear(GridStorage* storage) :
  storage(storage) {}

  /**
   * Destructor
   */
  ~OperationHierarchisationLinear() override {}

  void doHierarchisation(DataVector& node_values) override;
  void doDehierarchisation(DataVector& alpha) override;

 protected:
  /// Pointer to the grid's GridStorage object
  GridStorage* storage;
};

}  // namespace base
}  // namespace SGPP

#endif /* OPERATIONHIERARCHISATION_HPP */
