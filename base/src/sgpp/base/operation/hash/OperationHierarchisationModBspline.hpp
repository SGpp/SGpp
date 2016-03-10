// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONHIERARCHISATIONMODBSPLINE_HPP
#define OPERATIONHIERARCHISATIONMODBSPLINE_HPP

#include <sgpp/base/operation/hash/OperationHierarchisation.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/operation/hash/common/basis/BsplineModifiedBasis.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>


#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace base {

/**
 * Hierarchisation on sparse grid, mod bspline case
 *
 */
class OperationHierarchisationModBspline : public OperationHierarchisation {
 public:
  /**
   * Constructor
   *
   * @param storage the grid's GridStorage object
   * @param degree the bsplinenom's max. degree
   */
  explicit OperationHierarchisationModBspline(
    GridStorage& storage, size_t degree) : storage(storage), base(degree) {}

  /**
   * Destructor
   */
  ~OperationHierarchisationModBspline() override {}

  /**
   * Implements the hierarchisation on a sprase grid with mod bspline base functions
   *
   * @param node_values the functions values in the node base
   *
   */
  void doHierarchisation(DataVector& node_values) override;

  /**
   * Implements the dehierarchisation on a sprase grid with mod bspline base functions
   *
   * @param alpha the coefficients of the sparse grid's base functions
   *
   */
  void doDehierarchisation(DataVector& alpha) override;

 protected:
  /// Pointer to GridStorage object
  GridStorage& storage;
  /// Mod Bspline Basis object
  SBsplineModifiedBase base;
};

}  // namespace base
}  // namespace sgpp

#endif /* OPERATIONHIERARCHISATIONMODBSPLINE_HPP */
