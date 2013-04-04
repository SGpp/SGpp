/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Stefanie Schraufstetter (schraufs@in.tum.de)

#include "finance/basis/linear/boundary/algorithm_sweep/PhidPhiUpBBLinearBoundary.hpp"

namespace sg {
  namespace finance {



    PhidPhiUpBBLinearBoundary::PhidPhiUpBBLinearBoundary(sg::base::GridStorage* storage) : PhidPhiUpBBLinear(storage) {
    }

    PhidPhiUpBBLinearBoundary::~PhidPhiUpBBLinearBoundary() {
    }

    void PhidPhiUpBBLinearBoundary::operator()(sg::base::DataVector& source, sg::base::DataVector& result, grid_iterator& index, size_t dim) {
      // get boundary values
      double fl = 0.0;
      double fr = 0.0;

      // the following computations are independent from a bounding box
      if (!index.hint()) {
        index.top(dim);

        if (!this->storage->end(index.seq())) {
          rec(source, result, index, dim, fl, fr);
        }

        index.left_levelzero(dim);
      }

      size_t seq_left;
      size_t seq_right;

      // left boundary
      seq_left = index.seq();

      // right boundary
      index.right_levelzero(dim);
      seq_right = index.seq();

      // check boundary conditions
      if (this->boundingBox->hasDirichletBoundaryLeft(dim)) {
        result[seq_left] = 0.0; // source[seq_left];
      } else {
        // up
        //////////////////////////////////////
        result[seq_left] = fl;

        result[seq_left] += source[seq_right] * (-0.5);
      }

      if (this->boundingBox->hasDirichletBoundaryRight(dim)) {
        result[seq_right] = 0.0; //source[seq_right];
      } else {
        result[seq_right] = fr;
      }

      index.left_levelzero(dim);
    }

    // namespace detail

  } // namespace sg
}
