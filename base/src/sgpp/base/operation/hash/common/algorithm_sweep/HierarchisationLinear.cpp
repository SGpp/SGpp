// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/base/operation/hash/common/algorithm_sweep/HierarchisationLinear.hpp>

#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace base {



HierarchisationLinear::HierarchisationLinear(GridStorage& storage) : storage(
    storage) {
}

HierarchisationLinear::~HierarchisationLinear() {
}

void HierarchisationLinear::operator()(DataVector& source, DataVector& result,
                                       grid_iterator& index, size_t dim) {
  rec(source, result, index, dim, 0.0, 0.0);
}

void HierarchisationLinear::rec(DataVector& source, DataVector& result,
                                grid_iterator& index, size_t dim,
                                double fl, double fr) {
  // current position on the grid
  size_t seq = index.seq();
  // value in the middle, needed for recursive call and
  // calculation of the hierarchical surplus
  double fm = source[seq];

  // recursive calls for the right and left side of the current node
  if (index.hint() == false) {
    // descend left
    index.leftChild(dim);

    if (!storage.isInvalidSequenceNumber(index.seq())) {
      rec(source, result, index, dim, fl, fm);
    }

    // descend right
    index.stepRight(dim);

    if (!storage.isInvalidSequenceNumber(index.seq())) {
      rec(source, result, index, dim, fm, fr);
    }

    // ascend
    index.up(dim);
  }

  // hierarchisation
  result[seq] = fm - ((fl + fr) / 2.0);
}

}  // namespace base
}  // namespace sgpp
