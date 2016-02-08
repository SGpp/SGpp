// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef HIERARCHISATIONMODFUNDAMENTALSPLINE_HPP
#define HIERARCHISATIONMODFUNDAMENTALSPLINE_HPP

#include <sgpp/globaldef.hpp>
#include <sgpp/base/grid/type/ModFundamentalSplineGrid.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>

namespace SGPP {
namespace base {

/**
 * Functor for hierarchization with modified fundamental splines via
 * BreadthFirstSearch.
 */
class HierarchisationModFundamentalSpline {
 protected:
  /// grid iterator
  typedef GridStorage::grid_iterator grid_iterator;

 public:
  /**
   * Constructor.
   *
   * @param grid grid
   */
  HierarchisationModFundamentalSpline(ModFundamentalSplineGrid* grid);

  /**
   * Destructor.
   */
  virtual ~HierarchisationModFundamentalSpline();

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
  ModFundamentalSplineGrid* grid;
  /// grid storage
  GridStorage* storage;
};

}  // namespace base
}  // namespace SGPP

#endif /* HIERARCHISATIONMODFUNDAMENTALSPLINE_HPP */
