// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/basis/linearstretched/boundary/algorithm_sweep/SqXdPhidPhiDownBBLinearStretchedBoundary.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace finance {



SqXdPhidPhiDownBBLinearStretchedBoundary::SqXdPhidPhiDownBBLinearStretchedBoundary(
  SGPP::base::GridStorage* storage) : SqXdPhidPhiDownBBLinearStretched(storage) {
}

SqXdPhidPhiDownBBLinearStretchedBoundary::~SqXdPhidPhiDownBBLinearStretchedBoundary() {
}

void SqXdPhidPhiDownBBLinearStretchedBoundary::operator()(
  SGPP::base::DataVector& source, SGPP::base::DataVector& result,
  grid_iterator& index, size_t dim) {
  float_t q = this->stretching->getIntervalWidth(dim);
  float_t t = this->stretching->getIntervalOffset(dim);



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
  index.resetToLeftLevelZero(dim);
  seq_left = index.seq();
  left_boundary = source[seq_left];

  // right boundary
  index.resetToRightLevelZero(dim);
  seq_right = index.seq();
  right_boundary = source[seq_right];


  float_t bbFactor = ((q * q) + (3.0 * q * t) + (3.0 * t * t)) / (q);

  // check boundary conditions
  if (this->stretching->hasDirichletBoundaryLeft(dim)) {
    result[seq_left] = 0.0; //left_boundary;
  } else {
    result[seq_left] = (1.0 / 3.0) * left_boundary * bbFactor;
  }

  if (this->stretching->hasDirichletBoundaryRight(dim)) {
    result[seq_right] = 0.0; //right_boundary;
  } else {
    result[seq_right] = (1.0 / 3.0) * right_boundary * bbFactor;
    // down
    //////////////////////////////////////
    result[seq_right] -= (1.0 / 3.0) * left_boundary * bbFactor;
  }

  // move to root
  if (!index.hint()) {
    index.resetToLevelOne(dim);

    if (!this->storage->end(index.seq())) {
      rec(source, result, index, dim, left_boundary, right_boundary);
    }

    index.resetToLeftLevelZero(dim);
  }


}

// namespace detail

} // namespace SGPP
}