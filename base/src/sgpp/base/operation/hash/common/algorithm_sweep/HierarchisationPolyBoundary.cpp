// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/common/algorithm_sweep/HierarchisationPolyBoundary.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <cmath>

namespace sgpp {

namespace base {

HierarchisationPolyBoundary::HierarchisationPolyBoundary(GridStorage& storage,
                                                         SPolyBoundaryBase* base)
    : storage(storage), base(base) {}

HierarchisationPolyBoundary::~HierarchisationPolyBoundary() {}

void HierarchisationPolyBoundary::operator()(DataVector& source, DataVector& result,
                                             grid_iterator& index, size_t dim) {
  DataVector coeffs(index.getGridDepth(dim) + 2);

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

  // move to root
  if (!index.hint()) {
    index.resetToLevelOne(dim);

    if (!storage.isInvalidSequenceNumber(index.seq())) {
      rec(source, result, index, dim, coeffs);
    }

    index.resetToLeftLevelZero(dim);
  }
}

void HierarchisationPolyBoundary::rec(DataVector& source, DataVector& result, grid_iterator& index,
                                      size_t dim, DataVector& coeffs) {
  // current position on the grid
  size_t seq = index.seq();

  level_type cur_lev;
  index_type cur_ind;

  // get current level and index from grid
  index.get(dim, cur_lev, cur_ind);

  // hierarchisation
  double x = static_cast<double>(cur_ind) / static_cast<double>(1 << cur_lev);
  // v_i * 1 - sum_{j < i} v_j * \phi(x_i)
  result[seq] = source[seq] - base->evalHierToTop(cur_lev, cur_ind, coeffs, x);

  // recursive calls for the right and left side of the current node
  if (index.hint() == false) {
    // add the current result to the coefficient vector
    coeffs[cur_lev + 1] = result[seq];

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

    // remove the current result from the coefficient vector
    coeffs[cur_lev + 1] = 0.0;
  }
}

}  // namespace base
}  // namespace sgpp
