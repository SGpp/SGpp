/******************************************************************************
 * Copyright (C) 2011 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
// @author Dirk Pflueger (dirk.pflueger@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de), Sarpkan Selcuk (Sarpkan.Selcuk@mytum.de)

#include "pde/basis/linearstretched/boundary/algorithm_sweep/PhiPhiUpBBLinearStretchedBoundary.hpp"

namespace sg {
  namespace pde {



    PhiPhiUpBBLinearStretchedBoundary::PhiPhiUpBBLinearStretchedBoundary(sg::base::GridStorage* storage) : PhiPhiUpBBLinearStretched(storage) {
    }

    PhiPhiUpBBLinearStretchedBoundary::~PhiPhiUpBBLinearStretchedBoundary() {
    }

    void PhiPhiUpBBLinearStretchedBoundary::operator()(sg::base::DataVector& source, sg::base::DataVector& result, grid_iterator& index, size_t dim) {
      double q = this->stretching->getIntervalWidth(dim);
      //  double t = this->stretching->getIntervalOffset(dim);


      // get boundary values
      double fl = 0.0;
      double fr = 0.0;

      if (!index.hint()) {
        index.top(dim);

        if (!this->storage->end(index.seq())) {
          //This will be changed to rec
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
      //Left
      if (this->stretching->hasDirichletBoundaryLeft(dim)) {
        result[seq_left] = 0.0; // source[seq_left];
      } else {
        result[seq_left] = fl;
        result[seq_left] += (((1.0 / 6.0) * source[seq_right]) * q);
      }

      // Right
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
