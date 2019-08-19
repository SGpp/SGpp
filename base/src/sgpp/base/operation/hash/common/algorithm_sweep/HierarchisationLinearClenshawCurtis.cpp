// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/common/algorithm_sweep/HierarchisationLinearClenshawCurtis.hpp>

#include <sgpp/base/tools/ClenshawCurtisTable.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/globaldef.hpp>

#include <cmath>

namespace sgpp {
namespace base {

HierarchisationLinearClenshawCurtis::HierarchisationLinearClenshawCurtis(GridStorage& storage)
    : storage(storage), clenshawCurtisTable(ClenshawCurtisTable::getInstance()) {}

HierarchisationLinearClenshawCurtis::~HierarchisationLinearClenshawCurtis() {}

void HierarchisationLinearClenshawCurtis::operator()(DataVector& source, DataVector& result,
                                                     grid_iterator& index, size_t dim) {
  rec(source, result, index, dim, 0.0, 0.0, 1.0, 0.0);
}

void HierarchisationLinearClenshawCurtis::rec(DataVector& source, DataVector& result,
                                              grid_iterator& index, size_t dim, double xl,
                                              double fl, double xr, double fr) {
  // current position on the grid
  size_t seq = index.seq();

  level_type cur_lev;
  index_type cur_ind;

  // get current level and index from grid
  index.get(dim, cur_lev, cur_ind);

  // hierarchisation
  double xm = clenshawCurtisTable.getPoint(cur_lev, cur_ind);
  double fm = source[seq];
  result[seq] = fm - ((fr - fl) / (xr - xl) * (xm - xl) + fl);

  // recursive calls for the right and left side of the current node
  if (index.hint() == false) {
    // descend left
    index.leftChild(dim);

    if (!storage.isInvalidSequenceNumber(index.seq())) {
      rec(source, result, index, dim, xl, fl, xm, fm);
    }

    // descend right
    index.stepRight(dim);

    if (!storage.isInvalidSequenceNumber(index.seq())) {
      rec(source, result, index, dim, xm, fm, xr, fr);
    }

    // ascend
    index.up(dim);
  }
}

}  // namespace base
}  // namespace sgpp
