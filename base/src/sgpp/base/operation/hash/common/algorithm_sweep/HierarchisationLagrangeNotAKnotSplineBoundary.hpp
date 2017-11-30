// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef HIERARCHISATIONLAGRANGENOTAKNOTSPLINEBOUNDARY_HPP
#define HIERARCHISATIONLAGRANGENOTAKNOTSPLINEBOUNDARY_HPP

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

/**
 * Class that implements the hierarchisation on a Lagrange not-a-knot spline sparse grid with
 * boundaries. Therefore the ()operator has to be implement in order to use the sweep algorithm for
 * the grid traversal
 */
class HierarchisationLagrangeNotAKnotSplineBoundary {
 protected:
  typedef GridStorage::grid_iterator grid_iterator;
  GridStorage& storage;

 public:
  explicit HierarchisationLagrangeNotAKnotSplineBoundary(GridStorage& storage);
  virtual ~HierarchisationLagrangeNotAKnotSplineBoundary();

  void operator()(DataVector& source, DataVector& result,
                          grid_iterator& index, size_t dim);
};

}  // namespace base
}  // namespace sgpp

#endif /* HIERARCHISATIONLAGRANGENOTAKNOTSPLINEBOUNDARY_HPP */
