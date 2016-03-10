// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef HIERARCHISATIONFUNDAMENTALSPLINE_HPP
#define HIERARCHISATIONFUNDAMENTALSPLINE_HPP

#include <sgpp/globaldef.hpp>
#include <sgpp/base/grid/type/FundamentalSplineGrid.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>

namespace sgpp {
namespace base {

/**
 * Functor for hierarchization with fundamental splines via
 * BreadthFirstSearch.
 */
class HierarchisationFundamentalSpline {
 protected:
  /// grid iterator
  typedef GridStorage::grid_iterator grid_iterator;

 public:
  /**
   * Constructor.
   *
   * @param grid grid
   */
  explicit HierarchisationFundamentalSpline(FundamentalSplineGrid* grid);

  /**
   * Destructor.
   */
  virtual ~HierarchisationFundamentalSpline();

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
  FundamentalSplineGrid* grid;
  /// grid storage
  GridStorage& storage;
};

}  // namespace base
}  // namespace sgpp

#endif /* HIERARCHISATIONFUNDAMENTALSPLINE_HPP */
