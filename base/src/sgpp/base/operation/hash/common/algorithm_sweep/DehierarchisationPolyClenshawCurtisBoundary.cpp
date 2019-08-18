// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/common/algorithm_sweep/DehierarchisationPolyClenshawCurtisBoundary.hpp>

#include <sgpp/base/tools/ClenshawCurtisTable.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/globaldef.hpp>

#include <cmath>

namespace sgpp {
namespace base {

DehierarchisationPolyClenshawCurtisBoundary::DehierarchisationPolyClenshawCurtisBoundary(
    GridStorage& storage, SPolyClenshawCurtisBoundaryBase* base)
    : storage(storage), base(base), clenshawCurtisTable(ClenshawCurtisTable::getInstance()) {}

DehierarchisationPolyClenshawCurtisBoundary::~DehierarchisationPolyClenshawCurtisBoundary() {}

void DehierarchisationPolyClenshawCurtisBoundary::operator()(DataVector& source, DataVector& result,
                                                     grid_iterator& index, size_t dim) {
  DataVector coeffs(index.getGridDepth(dim) + 2);
  coeffs.setAll(0.0);

  // init coefficients at the boundary
  size_t seq;
  // left boundary
  index.resetToLeftLevelZero(dim);
  seq = index.seq();
  coeffs[0] = source[seq];
  // right boundary
  index.resetToRightLevelZero(dim);
  seq = index.seq();
  coeffs[1] = source[seq];

  rec(source, result, index, dim, coeffs);
}

void DehierarchisationPolyClenshawCurtisBoundary::rec(DataVector& source, DataVector& result,
                                              grid_iterator& index, size_t dim,
                                              DataVector& coeffs) {
  // current position on the grid
  size_t seq = index.seq();

  level_type cur_lev;
  index_type cur_ind;

  // get current level and index from grid
  index.get(dim, cur_lev, cur_ind);

  // Dehierarchisation
  // v_i * 1 + sum_{j < i} v_j * \phi(x_i)
  double x = clenshawCurtisTable.getPoint(cur_lev, cur_ind);
  result[seq] = source[seq] + base->evalHierToTop(cur_lev, cur_ind, coeffs, x);

  // recursive calls for the right and left side of the current node
  if (index.hint() == false) {
    // collect the hierarchical surpluses in the coefficient vector
    if (cur_lev == 0) {
      coeffs[cur_ind] = source[seq];
    } else {
      coeffs[cur_lev + 1] = source[seq];
    }

    // descend left
    index.leftChild(dim);

    if (!storage.isInvalidSequenceNumber(index.seq())) {
      rec(source, result, index, dim, coeffs);
    }

    // descend right
    index.stepRight(dim);

    if (!storage.isInvalidSequenceNumber(index.seq())) {
      rec(source, result, index, dim, coeffs);
    }

    // ascend
    index.up(dim);

    if (cur_lev == 0) {
      coeffs[cur_ind] = 0.0f;
    } else {
      coeffs[cur_lev + 1] = 0.0f;
    }
  }
}

}  // namespace base
}  // namespace sgpp
