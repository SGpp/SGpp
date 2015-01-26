/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Sam Maurus (MA thesis)

#include "finance/basis/linear/boundary/algorithm_sweep/XdPhidPhiDownBBLinearBoundary.hpp"
#include "base/exception/application_exception.hpp"

namespace sg {
  namespace finance {



    XdPhidPhiDownBBLinearBoundary::XdPhidPhiDownBBLinearBoundary(sg::base::GridStorage* storage) : XdPhidPhiDownBBLinear(storage) {
    }

    XdPhidPhiDownBBLinearBoundary::~XdPhidPhiDownBBLinearBoundary() {
    }

    void XdPhidPhiDownBBLinearBoundary::operator()(sg::base::DataVector& source, sg::base::DataVector& result, grid_iterator& index, size_t dim) {
      double q = this->boundingBox->getIntervalWidth(dim);
      double t = this->boundingBox->getIntervalOffset(dim);

      bool useBB = false;

      if (q != 1.0 || t != 0.0) {
        useBB = true;
      }

      // get boundary values
      double left_boundary;
      double right_boundary;
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

      if (useBB) {
        // check boundary conditions
        if (this->boundingBox->hasDirichletBoundaryLeft(dim)) {
          result[seq_left] = 0.0; //left_boundary;
        } else {
          throw new base::application_exception("XdPhidPhiDownBBLinearBoundary::operator : Not yet implemented for non-Dirichlet boundaries.");
        }

        if (this->boundingBox->hasDirichletBoundaryRight(dim)) {
          result[seq_right] = 0.0; //right_boundary;
        } else {
          throw new base::application_exception("XdPhidPhiDownBBLinearBoundary::operator : Not yet implemented for non-Dirichlet boundaries.");
        }

        // move to root
        if (!index.hint()) {
          index.top(dim);

          if (!this->storage->end(index.seq())) {
            recBB(source, result, index, dim, left_boundary, right_boundary, q, t);
          }

          index.left_levelzero(dim);
        }
      } else {
        // check boundary conditions
        if (this->boundingBox->hasDirichletBoundaryLeft(dim)) {
          result[seq_left] = 0.0; //left_boundary;
        } else {
          throw new base::application_exception("XdPhidPhiDownBBLinearBoundary::operator : Not yet implemented for non-Dirichlet boundaries.");
        }

        if (this->boundingBox->hasDirichletBoundaryRight(dim)) {
          result[seq_right] = 0.0; //right_boundary;
        } else {
          throw new base::application_exception("XdPhidPhiDownBBLinearBoundary::operator : Not yet implemented for non-Dirichlet boundaries.");
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
    }

    // namespace detail

  } // namespace sg
}
