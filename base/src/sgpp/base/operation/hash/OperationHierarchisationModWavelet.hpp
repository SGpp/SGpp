// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONHIERARCHISATIONMODWAVELET_HPP
#define OPERATIONHIERARCHISATIONMODWAVELET_HPP

#include <sgpp/base/operation/hash/OperationHierarchisation.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {

/**
 * Hierarchisation on sparse grid, mod wavelet case
 */
class OperationHierarchisationModWavelet : public OperationHierarchisation {
 public:
  /**
   * Constructor
   *
   * @param storage the grid's GridStorage object
   */
  explicit OperationHierarchisationModWavelet(GridStorage& storage) :
    storage(storage) {}

  /**
   * Destructor
   */
  ~OperationHierarchisationModWavelet() override {}

  /**
   * Implements the hierarchisation on a sprase grid with mod wavelets base functions
   *
   * @param node_values the functions values in the node base
   *
   */
  void doHierarchisation(DataVector& node_values) override;

  /**
   * Implements the dehierarchisation on a sprase grid with mod wavelets base functions
   *
   * @param alpha the coefficients of the sparse grid's base functions
   *
   */
  void doDehierarchisation(DataVector& alpha) override;

 protected:
  /// Pointer to GridStorage object
  GridStorage& storage;
};

}  // namespace base
}  // namespace SGPP

#endif /* OPERATIONHIERARCHISATIONMODWAVELET_HPP */
