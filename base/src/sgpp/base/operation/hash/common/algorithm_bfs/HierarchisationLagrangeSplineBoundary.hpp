// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef HIERARCHISATIONLAGRANGESPLINEBOUNDARY_HPP
#define HIERARCHISATIONLAGRANGESPLINEBOUNDARY_HPP

#include <sgpp/globaldef.hpp>
#include <sgpp/base/grid/type/LagrangeSplineBoundaryGrid.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>

namespace sgpp {
namespace base {

/**
 * Functor for hierarchization with fundamental splines via
 * BreadthFirstSearch.
 */
class HierarchisationLagrangeSplineBoundary {
 protected:
  /// grid iterator
  typedef GridStorage::grid_iterator grid_iterator;

 public:
  /**
   * Constructor.
   *
   * @param grid grid
   */
  explicit HierarchisationLagrangeSplineBoundary(LagrangeSplineBoundaryGrid* grid);

  /**
   * Destructor.
   */
  virtual ~HierarchisationLagrangeSplineBoundary();

  /**
   * Functor operator.
   * For each grid point, subtract the value of basis function at
   * the given iterator from the entry in result corresponding to the
   * grid point.
   *
   * @param[in]  source     node values
   * @param[out] result     result of the functor
   * @param      iterator   current grid point
   */
  void operator()(const DataVector& source,
                  DataVector& result,
                  const grid_iterator& iterator);

  /**
   * Functor operator.
   * For each grid point, subtract the value of basis function at
   * the given iterator from the row in result corresponding to the
   * grid point.
   *
   * @param[in]  source     node values
   * @param[out] result     result of the functor
   * @param      iterator   current grid point
   */
  void operator()(const DataMatrix& source,
                  DataMatrix& result,
                  const grid_iterator& iterator);

 protected:
  /// grid
  LagrangeSplineBoundaryGrid* grid;
  /// grid storage
  GridStorage& storage;
};

}  // namespace base
}  // namespace sgpp

#endif /* HIERARCHISATIONLAGRANGESPLINEBOUNDARY_HPP */
