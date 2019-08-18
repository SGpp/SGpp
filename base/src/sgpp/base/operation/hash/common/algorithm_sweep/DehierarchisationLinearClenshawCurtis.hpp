// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/base/tools/ClenshawCurtisTable.hpp>

#include <sgpp/globaldef.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearClenshawCurtisBasis.hpp>

namespace sgpp {
namespace base {

/**
 * Class that implements the dehierarchisation on a polynomial sparse grid. Therefore
 * the ()operator has to be implement in order to use the sweep algorithm for
 * the grid traversal
 */
class DehierarchisationLinearClenshawCurtis {
 protected:
  typedef GridStorage::grid_iterator grid_iterator;
  typedef level_t level_type;
  typedef index_t index_type;

  /// the grid object
  GridStorage& storage;

  /// clenshaw curtis points
  ClenshawCurtisTable& clenshawCurtisTable;

 public:
  /**
   * Constructor, must be bind to a grid
   *
   * @param storage the grid storage object of the the grid, on which the dehierarchisation should
   * be executed
   */
  explicit DehierarchisationLinearClenshawCurtis(GridStorage& storage);
  /**
   * Destructor
   */
  ~DehierarchisationLinearClenshawCurtis();

  /**
   * Implements operator() needed by the sweep class during the grid traversal. This function
   * is applied to the whole grid.
   *
   * @param source this DataVector holds the linear base coefficients of the sparse grid's
   * ansatz-functions
   * @param result this DataVector holds the node base coefficients of the function that should be
   * applied to the sparse grid
   * @param index a iterator object of the grid
   * @param dim current fixed dimension of the 'execution direction'
   */
  void operator()(DataVector& source, DataVector& result, grid_iterator& index, size_t dim);

 protected:
  /**
   * Recursive dehierarchisation algorithm, this algorithms works in-place -> source should be equal
   * to result
   *
   * @param source this DataVector holds the linear base coefficients of the sparse grid's
   * ansatz-functions
   * @param result this DataVector holds the node base coefficients of the function that should be
   * applied to the sparse grid
   * @param index a iterator object of the grid
   * @param dim current fixed dimension of the 'execution direction'
   * @param xl position of closest neighbor to the left
   * @param fl function value of closest neighbor to the left
   * @param xr position of closest neighbor to the right
   * @param fr function value of closest neighbor to the right
   */
  void rec(DataVector& source, DataVector& result, grid_iterator& index, size_t dim, double xl,
           double fl, double xr, double fr);
};

}  // namespace base
}  // namespace sgpp
