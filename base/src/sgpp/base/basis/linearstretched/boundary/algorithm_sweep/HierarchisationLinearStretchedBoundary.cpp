/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Sarpkan Selcuk (Sarpkan.Selcuk@mytum.de)

#include <sgpp/base/basis/linearstretched/boundary/algorithm_sweep/HierarchisationLinearStretchedBoundary.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {



    HierarchisationLinearStretchedBoundary::HierarchisationLinearStretchedBoundary(GridStorage* storage) : HierarchisationLinearStretched(storage) {
    }

    HierarchisationLinearStretchedBoundary::~HierarchisationLinearStretchedBoundary() {
    }

    void HierarchisationLinearStretchedBoundary::operator()(DataVector& source, DataVector& result, grid_iterator& index, size_t dim) {
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
