// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/base/operation/hash/common/algorithm_bfs/HierarchisationFundamentalNakSplineBoundary.hpp>

namespace sgpp {
namespace base {
HierarchisationFundamentalNakSplineBoundary::
HierarchisationFundamentalNakSplineBoundary(
  FundamentalNakSplineBoundaryGrid* grid) :
  grid(grid),
  storage(grid->getStorage()) {
}

HierarchisationFundamentalNakSplineBoundary::
~HierarchisationFundamentalNakSplineBoundary() {
}

void HierarchisationFundamentalNakSplineBoundary::operator()(
  const DataVector& source,
  DataVector& result,
  const grid_iterator& iterator) {
  const size_t n = storage.getSize();
  const size_t d = storage.getDimension();
  const size_t pointIndex = iterator.seq();

  SFundamentalNakSplineBase base(grid->getDegree());

  for (size_t q = 0; q < n; q++) {
    const GridPoint& point = storage[q];
    bool skipChild = false;

    if (q == pointIndex) {
      continue;
    }

    for (size_t t = 0; t < d; t++) {
      GridPoint::level_type l;
      GridPoint::level_type i;
      iterator.get(t, l, i);

      GridPoint::level_type k;
      GridPoint::level_type j;
      point.get(t, k, j);

      if ((k <= l) && ((k != l) || (i != j))) {
        skipChild = true;
        break;
      }
    }

    if (!skipChild) {
      double value = 1.0;

      for (size_t t = 0; t < d; t++) {
        GridPoint::level_type l;
        GridPoint::level_type i;
        iterator.get(t, l, i);

        const double val1d = base.eval(l, i, point.getStandardCoordinate(t));

        if (val1d == 0.0) {
          value = 0.0;
          break;
        }

        value *= val1d;
      }

      if (value != 0.0) {
        result[q] -= result[pointIndex] * value;
      }
    }
  }
}

void HierarchisationFundamentalNakSplineBoundary::operator()(
  const DataMatrix& source,
  DataMatrix& result,
  const grid_iterator& iterator) {
  const size_t n = storage.getSize();
  const size_t d = storage.getDimension();
  const size_t pointIndex = iterator.seq();

  SFundamentalNakSplineBase base(grid->getDegree());

  for (size_t q = 0; q < n; q++) {
    const GridPoint& point = storage[q];
    bool skipChild = false;

    if (q == pointIndex) {
      continue;
    }

    for (size_t t = 0; t < d; t++) {
      GridPoint::level_type l;
      GridPoint::level_type i;
      iterator.get(t, l, i);

      GridPoint::level_type k;
      GridPoint::level_type j;
      point.get(t, k, j);

      if ((k <= l) && ((k != l) || (i != j))) {
        skipChild = true;
        break;
      }
    }

    if (!skipChild) {
      double value = 1.0;

      for (size_t t = 0; t < d; t++) {
        GridPoint::level_type l;
        GridPoint::level_type i;
        iterator.get(t, l, i);

        const double val1d = base.eval(l, i, point.getStandardCoordinate(t));

        if (val1d == 0.0) {
          value = 0.0;
          break;
        }

        value *= val1d;
      }

      if (value != 0.0) {
        for (size_t j = 0; j < result.getNcols(); j++) {
          result.set(q, j, result.get(q, j) -
                     result.get(pointIndex, j) * value);
        }
      }
    }
  }
}

}  // namespace base
}  // namespace sgpp
