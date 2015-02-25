// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/basis/linear/boundary/algorithm_sweep/DPhiPhiDownBBLinearBoundary.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace finance {



    DPhiPhiDownBBLinearBoundary::DPhiPhiDownBBLinearBoundary(SGPP::base::GridStorage* storage) : DPhiPhiDownBBLinear(storage) {
    }

    DPhiPhiDownBBLinearBoundary::~DPhiPhiDownBBLinearBoundary() {
    }


    void DPhiPhiDownBBLinearBoundary::operator()(SGPP::base::DataVector& source, SGPP::base::DataVector& result, grid_iterator& index, size_t dim) {
      // get boundary values
      float_t left_boundary;
      float_t right_boundary;
      size_t seq_left;
      size_t seq_right;

      /*
       * Handle Level 0
       */
      // This handles the diagonal only
      //////////////////////////////////////
      // left boundary
      index.left_levelzero(dim);
      seq_left = index.seq();
      left_boundary = source[seq_left];

      // right boundary
      index.right_levelzero(dim);
      seq_right = index.seq();
      right_boundary = source[seq_right];


      // the following is independent from a bounding box

      // check boundary conditions
      if (this->boundingBox->hasDirichletBoundaryLeft(dim)) {
        result[seq_left] = 0.0; //left_boundary;
      } else {
        result[seq_left] = left_boundary * (-0.5);
      }

      if (this->boundingBox->hasDirichletBoundaryRight(dim)) {
        result[seq_right] = 0.0; //right_boundary;
      } else {
        result[seq_right] = right_boundary * (0.5);
        // down
        //////////////////////////////////////
        result[seq_right] += left_boundary * (-0.5);
      }

      // move to root
      if (!index.hint()) {
        index.top(dim);

        if (!this->storage->end(index.seq())) {
          rec(source, result, index, dim, left_boundary, right_boundary);
        }

        index.left_levelzero(dim);
      }
    }

    // namespace detail

  } // namespace SGPP
}