// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/common/algorithm_sweep/HierarchisationPolyClenshawCurtis.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/tools/ClenshawCurtisTable.hpp>

#include <sgpp/globaldef.hpp>

#include <cmath>

namespace sgpp {
namespace base {

HierarchisationPolyClenshawCurtis::HierarchisationPolyClenshawCurtis(GridStorage& storage,
                                                                     SPolyClenshawCurtisBase* base)
    : storage(storage), base(base), clenshawCurtisTable(ClenshawCurtisTable::getInstance()) {}

HierarchisationPolyClenshawCurtis::~HierarchisationPolyClenshawCurtis() {}

void HierarchisationPolyClenshawCurtis::operator()(DataVector& source, DataVector& result,
                                                   grid_iterator& index, size_t dim) {
  DataVector coeffs(index.getGridDepth(dim) + 1);

  rec(source, result, index, dim, coeffs);
}

void HierarchisationPolyClenshawCurtis::rec(DataVector& source, DataVector& result,
                                            grid_iterator& index, size_t dim, DataVector& coeffs) {
  // current position on the grid
  size_t seq = index.seq();

  level_type cur_lev;
  index_type cur_ind;

  // get current level and index from grid
  index.get(dim, cur_lev, cur_ind);

  // hierarchisation
  double x = clenshawCurtisTable.getPoint(cur_lev, cur_ind);
  result[seq] = source[seq] - base->evalHierToTop(cur_lev, cur_ind, coeffs, x);

  // recursive calls for the right and left side of the current node
  if (index.hint() == false) {
    coeffs[cur_lev] = result[seq];

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

    coeffs[cur_lev] = 0.0;
  }
}

}  // namespace base
}  // namespace sgpp
