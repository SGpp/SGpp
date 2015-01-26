/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include <sgpp/finance/basis/linear/boundary/algorithm_sweep/SqXdPhidPhiUpBBLinearBoundary.hpp>

namespace sg {
  namespace finance {



    SqXdPhidPhiUpBBLinearBoundary::SqXdPhidPhiUpBBLinearBoundary(sg::base::GridStorage* storage) : SqXdPhidPhiUpBBLinear(storage) {
    }

    SqXdPhidPhiUpBBLinearBoundary::~SqXdPhidPhiUpBBLinearBoundary() {
    }

    void SqXdPhidPhiUpBBLinearBoundary::operator()(sg::base::DataVector& source, sg::base::DataVector& result, grid_iterator& index, size_t dim) {
      double q = this->boundingBox->getIntervalWidth(dim);
      double t = this->boundingBox->getIntervalOffset(dim);

      bool useBB = false;

      if (q != 1.0 || t != 0.0) {
        useBB = true;
      }

      // get boundary values
      double fl = 0.0;
      double fr = 0.0;

      if (useBB) {
        if (!index.hint()) {
          index.top(dim);

          if (!this->storage->end(index.seq())) {
            recBB(source, result, index, dim, fl, fr, q, t);
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

        // up
        //////////////////////////////////////
        // check boundary conditions
        if (this->boundingBox->hasDirichletBoundaryLeft(dim)) {
          result[seq_left] = 0.0; // source[seq_left];
        } else {
          result[seq_left] = fl;
          double bbFactor = ((q * q) + (3.0 * q * t) + (3.0 * t * t)) / (q);
          result[seq_left] -= (1.0 / 3.0) * source[seq_right] * bbFactor;
        }

        if (this->boundingBox->hasDirichletBoundaryRight(dim)) {
          result[seq_right] = 0.0; //source[seq_right];
        } else {
          result[seq_right] = fr;
        }

        index.left_levelzero(dim);
      } else {
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

        // up
        //////////////////////////////////////
        // check boundary conditions
        if (this->boundingBox->hasDirichletBoundaryLeft(dim)) {
          result[seq_left] = 0.0; // source[seq_left];
        } else {
          result[seq_left] = fl;

          result[seq_left] -= (1.0 / 3.0) * source[seq_right];
        }

        if (this->boundingBox->hasDirichletBoundaryRight(dim)) {
          result[seq_right] = 0.0; //source[seq_right];
        } else {
          result[seq_right] = fr;
        }

        index.left_levelzero(dim);
      }
    }

    // namespace detail

  } // namespace sg
}
