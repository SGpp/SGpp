// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#include <sgpp/finance/basis/linearstretched/boundary/algorithm_sweep/PhidPhiUpBBLinearStretchedBoundary.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace finance {



    PhidPhiUpBBLinearStretchedBoundary::PhidPhiUpBBLinearStretchedBoundary(SGPP::base::GridStorage* storage) : PhidPhiUpBBLinearStretched(storage) {
    }

    PhidPhiUpBBLinearStretchedBoundary::~PhidPhiUpBBLinearStretchedBoundary() {
    }

    void PhidPhiUpBBLinearStretchedBoundary::operator()(SGPP::base::DataVector& source, SGPP::base::DataVector& result, grid_iterator& index, size_t dim) {
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
      if (this->stretching->hasDirichletBoundaryLeft(dim)) {
        result[seq_left] = 0.0; // source[seq_left];
      } else {
        // up
        //////////////////////////////////////
        result[seq_left] = fl;

        result[seq_left] += source[seq_right] * (-0.5);
      }

      if (this->stretching->hasDirichletBoundaryRight(dim)) {
        result[seq_right] = 0.0; //source[seq_right];
      } else {
        result[seq_right] = fr;
      }

      index.left_levelzero(dim);
    }

    // namespace detail

  } // namespace SGPP
}