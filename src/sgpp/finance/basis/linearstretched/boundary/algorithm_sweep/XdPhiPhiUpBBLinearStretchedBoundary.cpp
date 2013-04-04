/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Sarpkan Selcuk (Sarpkan.Selcuk@mytum.de)

#include "finance/basis/linearstretched/boundary/algorithm_sweep/XdPhiPhiUpBBLinearStretchedBoundary.hpp"

namespace sg {
  namespace finance {



    XdPhiPhiUpBBLinearStretchedBoundary::XdPhiPhiUpBBLinearStretchedBoundary(sg::base::GridStorage* storage) : XdPhiPhiUpBBLinearStretched(storage) {
    }

    XdPhiPhiUpBBLinearStretchedBoundary::~XdPhiPhiUpBBLinearStretchedBoundary() {
    }

    void XdPhiPhiUpBBLinearStretchedBoundary::operator()(sg::base::DataVector& source, sg::base::DataVector& result, grid_iterator& index, size_t dim) {
      double q = this->stretching->getIntervalWidth(dim);
      double t = this->stretching->getIntervalOffset(dim);



      // get boundary values
      double fl = 0.0;
      double fr = 0.0;


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

        result[seq_left] += (source[seq_right] * (((1.0 / 6.0) * q) + (0.5 * t)));
      }

      if (this->stretching->hasDirichletBoundaryRight(dim)) {
        result[seq_right] = 0.0; //source[seq_right];
      } else {
        result[seq_right] = fr;
      }

      index.left_levelzero(dim);



    }
    // namespace detail

  } // namespace sg
}
