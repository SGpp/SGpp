// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONHIERARCHISATIONMODLINEAR_HPP
#define OPERATIONHIERARCHISATIONMODLINEAR_HPP

#include <sgpp/base/operation/hash/OperationHierarchisation.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {

/**
 * Hierarchisation on sparse grid, mod linear case
 *
 */
class OperationHierarchisationModLinear : public OperationHierarchisation {
 public:
  /**
   * Constructor
   *
   * @param storage the grid's GridStorage object
   */
  OperationHierarchisationModLinear(GridStorage* storage) : storage(storage) {}

  /**
   * Destructor
   */
  virtual ~OperationHierarchisationModLinear() override {}

  virtual void doHierarchisation(DataVector& node_values) override;
  virtual void doDehierarchisation(DataVector& alpha) override;

 protected:
  /// Pointer to GridStorage object
  GridStorage* storage;
};

}  // namespace base
}  // namespace SGPP

#endif /* OPERATIONHIERARCHISATIONMODLINEAR_HPP */