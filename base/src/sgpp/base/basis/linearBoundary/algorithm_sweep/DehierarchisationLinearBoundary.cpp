// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#include <sgpp/base/basis/linearBoundary/algorithm_sweep/DehierarchisationLinearBoundary.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {



    DehierarchisationLinearBoundary::DehierarchisationLinearBoundary(GridStorage* storage) : DehierarchisationLinear(storage) {
    }

    DehierarchisationLinearBoundary::~DehierarchisationLinearBoundary() {
    }

    void DehierarchisationLinearBoundary::operator()(DataVector& source, DataVector& result, grid_iterator& index, size_t dim) {
      double left_boundary;
      double right_boundary;
      size_t seq;

      // left boundary
      index.left_levelzero(dim);
      seq = index.seq();
      left_boundary = source[seq];
      // right boundary
      index.right_levelzero(dim);
      seq = index.seq();
      right_boundary = source[seq];

      // move to root
      if (!index.hint()) {
        index.top(dim);

        if (!storage->end(index.seq())) {
          rec(source, result, index, dim, left_boundary, right_boundary);
        }

        index.left_levelzero(dim);
      }
    }

    // namespace detail

  } // namespace SGPP
}