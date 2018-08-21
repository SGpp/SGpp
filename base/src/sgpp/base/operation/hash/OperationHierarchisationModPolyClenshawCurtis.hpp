// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/operation/hash/OperationHierarchisation.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/operation/hash/common/basis/PolyModifiedClenshawCurtisBasis.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

/**
 * Hierarchisation on sparse grid, poly case
 */
class OperationHierarchisationModPolyClenshawCurtis : public OperationHierarchisation {
 public:
  /**
   * Constructor
   *
   * @param storage the grid's GridStorage object
   * @param degree the polynom's max. degree
   */
  OperationHierarchisationModPolyClenshawCurtis(GridStorage& storage, size_t degree)
      : storage(storage), base(degree) {}

  /**
   * Destructor
   */
  ~OperationHierarchisationModPolyClenshawCurtis() override {}

  /**
   * Implements the hierarchisation on a sprase grid with poly base functions
   *
   * @param node_values the functions values in the node base
   *
   */
  void doHierarchisation(DataVector& node_values) override;

  /**
   * Implements the dehierarchisation on a sprase grid with poly base functions
   *
   * @param alpha the coefficients of the sparse grid's base functions
   *
   */
  void doDehierarchisation(DataVector& alpha) override;

 protected:
  /// Pointer to GridStorage object
  GridStorage& storage;
  /// ModPoly Basis object
  SPolyModifiedClenshawCurtisBase base;
};

}  // namespace base
}  // namespace sgpp
