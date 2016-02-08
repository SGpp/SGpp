// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/common/algorithm_sweep/DehierarchisationLinearStretchedBoundary.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {



DehierarchisationLinearStretchedBoundary::DehierarchisationLinearStretchedBoundary(
  GridStorage* storage) : DehierarchisationLinearStretched(storage) {
}

DehierarchisationLinearStretchedBoundary::~DehierarchisationLinearStretchedBoundary() {
}

void DehierarchisationLinearStretchedBoundary::operator()(DataVector& source,
    DataVector& result, grid_iterator& index, size_t dim) {
  float_t left_boundary;
  float_t right_boundary;
  size_t seq;

  // left boundary
  index.resetToLeftLevelZero(dim);
  seq = index.seq();
  left_boundary = source[seq];
  // right boundary
  index.resetToRightLevelZero(dim);
  seq = index.seq();
  right_boundary = source[seq];

  // move to root
  if (!index.hint()) {
    index.resetToLevelOne(dim);

    if (!storage->end(index.seq())) {
      rec(source, result, index, dim, left_boundary, right_boundary);
    }

    index.resetToLeftLevelZero(dim);
  }
}

// namespace detail

} // namespace SGPP
}