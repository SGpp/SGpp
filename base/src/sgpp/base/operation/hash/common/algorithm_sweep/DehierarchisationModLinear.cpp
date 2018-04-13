// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/base/operation/hash/common/algorithm_sweep/DehierarchisationModLinear.hpp>

#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace base {



DehierarchisationModLinear::DehierarchisationModLinear(GridStorage& storage) :
  storage(storage) {
}

DehierarchisationModLinear::~DehierarchisationModLinear() {
}

void DehierarchisationModLinear::operator()(DataVector& source,
    DataVector& result, grid_iterator& index, size_t dim) {
  rec(source, result, index, dim, 0.0, 0.0);
}

void DehierarchisationModLinear::rec(DataVector& source, DataVector& result,
                                     grid_iterator& index, size_t dim,
                                     double fl, double fr) {
  // current position on the grid
  size_t seq = index.seq();
  // value in the middle, needed for recursive call and
  // calculation of the hierarchical surplus
  double fm = source[seq];

  // dehierarchisation
  fm += ((fl + fr) / 2.0);
  result[seq] = fm;

  level_t l;
  index_t i;

  index.get(dim, l, i);

  // recursive calls for the right and left side of the current node
  if (index.hint() == false) {
    double fltemp = fl;
    double frtemp = fr;

    // When we descend the hierarchical basis
    // we have to modify the boundary values
    // in case the index is 1 or (2^l)-1 or we are on the first level
    // level 1, constant function
    if (l == 1) {
      // constant function
      fltemp = fm;
      frtemp = fm;
    } else if (i == 1) {  // left boundary
      double ftemp;
      ftemp = fr - fm;
      fltemp = fm - ftemp;
    } else if (static_cast<int>(i) == static_cast<int>((1 << l) - 1)) {
      // right boundary
      double ftemp;
      ftemp = fl - fm;
      frtemp = fm - ftemp;
    } else {  // inner functions
    }

    // descend left
    index.leftChild(dim);

    if (!storage.isInvalidSequenceNumber(index.seq())) {
      rec(source, result, index, dim, fltemp, fm);
    }

    // descend right
    index.stepRight(dim);

    if (!storage.isInvalidSequenceNumber(index.seq())) {
      rec(source, result, index, dim, fm, frtemp);
    }

    // ascend
    index.up(dim);
  }
}

}  // namespace base
}  // namespace sgpp
