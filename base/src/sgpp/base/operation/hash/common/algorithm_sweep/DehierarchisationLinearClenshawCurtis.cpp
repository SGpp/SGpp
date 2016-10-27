// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/tools/ClenshawCurtisTable.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/globaldef.hpp>

#include <cmath>

#include "DehierarchisationLinearClenshawCurtis.hpp"

namespace sgpp {
namespace base {

DehierarchisationLinearClenshawCurtis::DehierarchisationLinearClenshawCurtis(GridStorage& storage)
    : storage(storage), clenshawCurtisTable(ClenshawCurtisTable::getInstance()) {}

DehierarchisationLinearClenshawCurtis::~DehierarchisationLinearClenshawCurtis() {}

void DehierarchisationLinearClenshawCurtis::operator()(DataVector& source, DataVector& result,
                                                       grid_iterator& index, size_t dim) {
  rec(source, result, index, dim, 0.0, 0.0, 1.0, 0.0);
}

void DehierarchisationLinearClenshawCurtis::rec(DataVector& source, DataVector& result,
                                                grid_iterator& index, size_t dim, double xl,
                                                double fl, double xr, double fr) {
  // current position on the grid
  size_t seq = index.seq();

  level_type cur_lev;
  index_type cur_ind;

  // get current level and index from grid
  index.get(dim, cur_lev, cur_ind);

  // Dehierarchisation
  // v_i * 1 + sum_{j < i} v_j * \phi(x_i)
  double xm = clenshawCurtisTable.getPoint(cur_lev, cur_ind);
  double fm = source[seq];
  result[seq] = fm + ((fr - fl) / (xr - xl) * (xm - xl) + fl);

  // recursive calls for the right and left side of the current node
  if (index.hint() == false) {
    // descend left
    index.leftChild(dim);

    if (!storage.isValidSequenceNumber(index.seq())) {
      rec(source, result, index, dim, xl, fl, xm, fm);
    }

    // descend right
    index.stepRight(dim);

    if (!storage.isValidSequenceNumber(index.seq())) {
      rec(source, result, index, dim, xm, fm, xr, fr);
    }

    // ascend
    index.up(dim);
  }
}

}  // namespace base
}  // namespace sgpp
