// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/operation/hash/OperationHierarchisation.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/grid/Grid.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

/**
 * This class implements the hierarchisation and dehierarchisation on sparse grids
 * which do not have a grid point at the boundary in each direction for every inner
 * point (boundaryLevel > 1)
 */
class OperationArbitraryBoundaryHierarchisation : public OperationHierarchisation {
 public:
  /**
   * Constructor
   */
  OperationArbitraryBoundaryHierarchisation(Grid& grid);

  /**
   * Destructor
   */
  virtual ~OperationArbitraryBoundaryHierarchisation();

  /**
   * Implements the hierarchisation on a sparse grid
   *
   * @param node_values the function's values in the nodal basis
   */
  void doHierarchisation(DataVector& nodal_values);

  /**
   * Implements the dehierarchisation on a sparse grid
   *
   * @param alpha the coefficients of the sparse grid's basis functions
   */
  void doDehierarchisation(DataVector& alpha);

 private:
  Grid& grid;
  std::unique_ptr<Grid> boundaryGrid;
  std::unique_ptr<Grid> innerGrid;

  OperationMultipleEval* createOperationMultipleEval(Grid& grid, DataMatrix& coordinates);
};

}  // namespace base
}  // namespace sgpp
