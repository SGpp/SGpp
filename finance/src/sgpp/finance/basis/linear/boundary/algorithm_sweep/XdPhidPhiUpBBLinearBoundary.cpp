/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Sam Maurus (MA thesis)

#include <sgpp/finance/basis/linear/boundary/algorithm_sweep/XdPhidPhiUpBBLinearBoundary.hpp>
#include <sgpp/base/exception/application_exception.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace finance {



    XdPhidPhiUpBBLinearBoundary::XdPhidPhiUpBBLinearBoundary(SGPP::base::GridStorage* storage) : XdPhidPhiUpBBLinear(storage) {
    }

    XdPhidPhiUpBBLinearBoundary::~XdPhidPhiUpBBLinearBoundary() {
    }

    void XdPhidPhiUpBBLinearBoundary::operator()(SGPP::base::DataVector& source, SGPP::base::DataVector& result, grid_iterator& index, size_t dim) {
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

        // check boundary conditions
        if (this->boundingBox->hasDirichletBoundaryLeft(dim)) {
          result[seq_left] = 0.0; // source[seq_left];
        } else {
          throw new base::application_exception("XdPhidPhiUpBBLinearBoundary::operator : Not yet implemented for non-Dirichlet boundaries.");
        }

        if (this->boundingBox->hasDirichletBoundaryRight(dim)) {
          result[seq_right] = 0.0; //source[seq_right];
        } else {
          throw new base::application_exception("XdPhidPhiUpBBLinearBoundary::operator : Not yet implemented for non-Dirichlet boundaries.");
        }

        index.left_levelzero(dim);
      } else {
        if (!index.hint()) {
          index.top(dim);

          if (!storage->end(index.seq())) {
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
          throw new base::application_exception("XdPhidPhiUpBBLinearBoundary::operator : Not yet implemented for non-Dirichlet boundaries.");
        }

        if (this->boundingBox->hasDirichletBoundaryRight(dim)) {
          result[seq_right] = 0.0; //source[seq_right];
        } else {
          throw new base::application_exception("XdPhidPhiUpBBLinearBoundary::operator : Not yet implemented for non-Dirichlet boundaries.");
        }

        index.left_levelzero(dim);
      }
    }

    // namespace detail

  } // namespace SGPP
}
